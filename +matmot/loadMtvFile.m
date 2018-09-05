function [data, meta] = loadMtvFile(filename, varargin)
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

import matmot.FormatSpec.*

inp = inputParser();
inp.addParameter('fixTimestamps', true);
inp.addParameter('fixTimestampsIncrement', 1e3);
inp.addParameter('fixTimestampsThresholdFrames', 120);
inp.parse(varargin{:});
P = inp.Results;

% Open the binary data file
[fid, msg] = fopen(filename, 'r');
if fid == -1
    error('Could not open file %s: message "%s"', filename, msg);
end

% Read off header (pre-v0.2 only)
[meta, ver, offset] = readMeta(filename);
fread(fid, offset);

nRb = meta.n_rigid_bodies;
nMrk = meta.n_markers;
nBytesFrame = bytesPerFrame(nRb, nMrk);
allBytes = fread(fid, [nBytesFrame, meta.n_frames], '*uint8');
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECODE FIELDS
% The "basic" fields occur once in each record. The remainder of fields
% correspond to rigid-body or marker data, and therefore are repeated 
% according to the number of rigid bodies and markers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get names, encodings and repetition number for all fields
allFields = [
    basicFields()
    rbFields()
    markerFields()
    ];

fieldNReps = [ 
    ones(size(basicFields())) 
    ones(size(rbFields()))*nRb
    ones(size(markerFields()))*nMrk
    ];


% Decode bytes for each field
startByteIdx = 1;

for f = 1:numel(allFields)
    field = allFields(f);
    nrep = fieldNReps(f);
    nb = field.n_bytes;
    inds = startByteIdx + (0 : (nrep-1))*nb + (0 : (nb-1))';
    encoding = convertEncoding(field.encoding);
    bytes = allBytes(inds(:), :);
    tmp = typecast(bytes(:), encoding);
    data.(field.name) = reshape(tmp, nrep, meta.n_frames)';
    startByteIdx = startByteIdx + nb*nrep;
end

if P.fixTimestamps
    data.frameTimestamp = fixTimestamps(data.frameTimestamp, P);
end

end

function t = fixTimestamps(t0, P)
% Ensures timestamps are monotonically ascending

t = t0;

increment = P.fixTimestampsIncrement;
thresh = P.fixTimestampsThresholdFrames;

dt = diff(t0);
iNeg = find(dt<0);
iNegBef = max(1, iNeg-thresh+1);
iNegAft = min(numel(t0), iNeg+thresh);

iBad = [];
for n = 1:numel(iNeg)
    indsBef = iNegBef(n):iNeg(n);
    indsAft = (iNeg(n)+1 : iNegAft(n));
    if min(t0(indsBef)) > max(t0(indsAft))
        iBad(end+1, 1) = iNeg(n);
    end
end

nBad = numel(iBad);
nFrames = numel(t0);

if nBad > 0
    fprintf('Timestamps are not monotonically ascending. Fixing %u negative jumps.\n', nBad);
    for i = 1:nBad
        inds = (iBad(i)+1 : nFrames)';
        tShift = t(inds(1)-1) - t(inds(1)) + increment;
        t(inds) = t(inds) + tShift;
    end
end

end