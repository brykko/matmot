function [data] = loadMtvFile(fileName)
%LOADMTVFILE read contents of binary .mtv file from Streamer
%
% DATA = LOADMTVFILE(FILENAME) reads the contents of the .mtv file
% specified by FILENAME, returning structure DATA. Fields of DATA are as
% follows:
%
%   -----------------------------------------------------------------------
%   FRAME DATA
%   These fields contain data about the motion capture frame. Each field's
%   data is a column vector with one element per frame.
%   
%   frameIdx (int32) - Motive frame ID. Frames are counted from when Motive
%   was opened.
%
%   frameTimestamp (double) - time of Motive frame, relative to when Motive
%   was opened.
%
%   frameLatency (single) - latency of Motive frame (not currently used).
%
%   -----------------------------------------------------------------------
%   RIGID BODY DATA
%   These fields contain data about the rigid body objects configured in 
%   Motive. Each field contains a matrix, with frames in rows and rigid
%   bodies in columns. The number of columns will equal the value of the
%   Streamer parameter "nRigidBodies" that was used for acquiring the data.
%
%   rbx, rby, rbz (single) - position coordinates of rigid body.
%
%   rbqx, rbqy, rbqz, rbqw (single) - quaternion rotations of rigid body.
%
%   posError (single) mean rigid body position error.
%
%   posTracked (uint8) rigid body tracking success (1 = tracked, 0 =
%   untracked)
%
%   -----------------------------------------------------------------------
%   MARKER DATA
%   These fields contain data about labelled markers. Each field contains a
%   matrix, with frames in rows and markers in columns. The number of 
%   columns will equal the value of the Streamer parameter "nMarkers" used 
%   for acquiring the data.
%
%   mx, my, mz (single) - coordinates of labelled markers.
%
%   msz (single) - size of labelled markers.
%
%   mres (single) - residuals of labelled markers.

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

header = fread(fid, [1, FormatSpec.HEADER_LENGTH], '*char');
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
