function [data] = loadMtvFile(fileName)
%LOADMTVFILE read contents of binary .mtv file from Streamer
%
% DATA = LOADMTVFILE(FILENAME) reads the contents of the .DAT file
% specified by FILENAME, return structure DATA. Fields of DATA are as
% follows:
%
%   frameIdx (int32) - Motive frame ID. Frames are counted from when Motive
%   was opened.
%
%   frameTimestamp (double) - time of Motive frame, relative to when Motive
%   was opened.
%
%   frameLatency (single) - latency of Motive frame.
%
%   pos (single) - [x, y, z] coordinates of rigid body
%
%   rot (single) - [qx, qy, qz, qw] quaternion rotations of rigid body
%
%   posError (single) mean rigid body position error
%
%   posTracked (uint8) rigid body tracking success (1 = tracked, 0 =
%   untracked)
%
%   mx, my, mz (single) - coordinates of labelled markers. These may
%   include markers included in identified rigid bodies. Each field
%   contains an nframes-by-nmarkers matrix, where nmarkers is the maximum
%   number of markers recorded, as set by parameter 'nMarkers' in
%   Streamer, represented by field "n_markers" in the .mtv file
%   header.
%
%   msz (single) - size of labelled markers. Format is the same as for mx,
%   my, mz.
%
%   mres (single) - residuals of labelled markers. Format is the same as
%   for mx, my, mz.

import matmot.Streamer
import matmot.FormatSpec

% Read all data from file
[pth, fn, ext] = fileparts(fileName);
if isempty(pth), pth = pwd(); end
fpth = fullfile(pth, [fn ext]);
assert(exist(fpth, 'file') > 0, 'File "%s" does not exist', fpth);
[fid, msg] = fopen(fpth, 'r');
if fid == -1
    error('Could not open file %s: message "%s"', fpth, msg);
end

header = fread(fid, [1, Streamer.HEADER_LENGTH], '*char');
fileVersion = getVersion(header);

if fileVersion(2) > 0
    nRbs = getHeaderField(header, 'n_rigid_bodies');
else
    nRbs = 1;
end
nMarkers = getHeaderField(header, 'n_markers');

nBytesBasic = FormatSpec.bytesBasic();
nBytesRb = nRbs * FormatSpec.bytesPerRb();
nBytesFrame = FormatSpec.bytesPerFrame(nRbs, nMarkers);
allBytes = fread(fid, [nBytesFrame, inf], '*uint8');
fclose(fid);

nFrames = size(allBytes, 2);
basicFields = FormatSpec.basicFields();
rbFields = FormatSpec.rbFields();

% Cycle through rigid body fields and interpret the relevant bytes for field
% appropriately
for f = 1:numel(basicFields)
    field = basicFields(f);
    inds = field.byte_inds + (0 : (field.n_bytes-1));
    bytes = allBytes(inds, :);
    tmp = typecast(bytes(:), field.encoding);
    data.(field.name) = reshape(tmp, [], nFrames)';
end

% Read the rigid bodies
indField = nBytesBasic + 1;
for f = 1:numel(rbFields)
    field = rbFields(f);
    nBytes = field.n_bytes;
    for r = 1:nRbs
        inds = indField + (r-1)*nBytes + (0 : (nBytes-1));
        bytes = allBytes(inds, :);
        tmp = typecast(bytes(:), field.encoding);
        data.(field.name)(:, r) = reshape(tmp, [], nFrames)';
    end
    indField = indField + field.n_bytes*nRbs;
end

% Read the markers
markerFields = FormatSpec.markerFields();
indField = nBytesBasic + nBytesRb + 1;
for f = 1:numel(markerFields)
    field = markerFields(f);
    encodingConv = FormatSpec.convertEncoding(field.encoding);
    nBytes = field.n_bytes;
    for m = 1:nMarkers
        inds = indField + (m-1)*nBytes + (0 : (nBytes-1));
        bytes = allBytes(inds, :);
        tmp = typecast(bytes(:), encodingConv);
        data.(field.name)(:, m) = reshape(tmp, [], nFrames)';
    end
    indField = indField + field.n_bytes*nMarkers;
end

fprintf('Done.\n');

end

function val = getHeaderField(header, field, fmt)

lines = strsplit(header, '\r\n');
expr = [field '=(\w+).*?'];
[mat, tok] = regexp(lines, expr, 'match', 'tokens');
idx = find(~cellfun(@isempty, mat));
assert(~isempty(idx), ...
    'importMtv:headerFieldNotFound', ...
    'Error parsing file header: could not find string "%s"', field)
val = tok{idx}{1}{1};
try
    % Convert to numeric if possible, otherwise leave as char
    val = str2num(val);
end

end

function v = getVersion(header)
lines = strsplit(header, '\r\n');
strs = strsplit(lines{1});
verStr = strs{end};
verElStr = strsplit(verStr, '.');
assert(numel(verElStr)==3, 'matmot:loadMtvFile:noFileVersion', ...
    'Failed to interpret version of .mtv file.');
for i = 1:3
    v(i) = str2num(verElStr{i});
end
end