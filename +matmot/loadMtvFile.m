function [data, meta] = loadMtvFile(filename)
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

% Determine file format version
% v0.1 has a single data file with 16 kb header
% v0.2+ has a headerless data file with separate metadata text file
txt = FormatSpec.readMetaText(filename);
meta = FormatSpec.parseMetaText(txt);
ver = meta.matmot_version;

if ver(2) >= 1
    nRbs = meta.n_rigid_bodies;
else
    nRbs = 1;
    meta.n_rigid_bodies = nRbs;
end

% Open the binary data file
[fid, msg] = fopen(filename, 'r');
if fid == -1
    error('Could not open file %s: message "%s"', filename, msg);
end

if ver(2) < 2
    fread(fid, FormatSpec.HEADER_LENGTH);
end

nBytesBasic = FormatSpec.bytesBasic();
nBytesRb = nRbs * FormatSpec.bytesPerRb();
nBytesFrame = FormatSpec.bytesPerFrame(nRbs, meta.n_markers);
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
    for m = 1:meta.n_markers
        inds = indField + (m-1)*nBytes + (0 : (nBytes-1));
        bytes = allBytes(inds, :);
        tmp = typecast(bytes(:), encodingConv);
        data.(field.name)(:, m) = reshape(tmp, [], nFrames)';
    end
    indField = indField + field.n_bytes*meta.n_markers;
end

end
