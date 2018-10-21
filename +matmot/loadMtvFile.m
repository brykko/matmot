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
inp.addParameter('frameCleanupMethod', 'sort');
inp.addParameter('fixTimestamps', []);
inp.addParameter('fixTimestampsIncrement', 1e3);
inp.addParameter('fixTimestampsThresholdFrames', 120);
inp.parse(varargin{:});
P = inp.Results;

if ~isempty(P.fixTimestamps)
    warning('matmot:loadMtvFile:fixTimestampsDeprecated', ...
        'Parameter "fixTimestamps" is deprecated, please use the new alternative "frameCleanupMethod" instead.');
    P.frameCleanupMethod = 'fixTimestamps';
end
    

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

% Check for bad frames
nt = numel(data.frameTimestamp);
[iNeg, iOverlap, iMatch] = checkTimestamps(data.frameTimestamp);
nNeg = numel(iNeg);
nOverlap = numel(iOverlap);
percentOverlap = nOverlap ./ (nt-nOverlap) .* 100;
nMatch = numel(find(~isnan(iMatch)));
percentDupTime = nMatch ./ nOverlap .* 100;
if nOverlap > 0
    warning( ...
        'matmot:loadMtvFile:overlappingTimestamps', ...
        'The timestamps in this file are not monotonically ascending (%u negative jumps, %.2f%% overlapping frames, of which %.2f%% have matching (duplicated) timestamps)', ...
        nNeg, percentOverlap, percentDupTime);
end

if strcmpi(P.frameCleanupMethod, 'fixTimestamps')
    % OLD METHOD: out-of-order timestamp frames are interpreted as having
    % erroneous timestamps, but being otherwise OK. All frames are kept (in
    % their original order), and timestamp values are altered to remove
    % overlap.
    data = fixTimestamps(data, P);
elseif any(strcmpi(P.frameCleanupMethod, {'sort', 'discard'}))
    % NEW METHOD: out-of-order timestamps are interpreted as representing
    % out-of-order frames which may or may not be duplicates of previous
    % frames in the file. The "discard" method eliminates all frames
    % occurring in overlapping time regions (most conservative), while the
    % "sort" method eliminates only frames whose timestamp exactly matches
    % a previous frame in the file.
    data = fixOverlappingFrames(data, iOverlap, iMatch, P.frameCleanupMethod);
elseif ~strcmpi(P.frameCleanupMethod, 'none')
    error('matmot:loadMtvFile:invalidFrameCleanupMethod', ...
        'Valid options for parameter "frameCleanupMethod" are "sort", "discard", "fixTimestamps" and "none"');
end

end

function [iNeg, iOverlap, iMatch] = checkTimestamps(t)

    dt = diff(t);
    iNeg = find(dt<0);
    
    lastValidTime = -inf;
    iOverlap = [];
    iMatch = [];
    nOverlap = 0;
    
    for n = 2:numel(t)
        
        if dt(n-1) < 0
            inOverlapZone = true;
        else
            inOverlapZone = t(n) <= lastValidTime;
        end
        
        if inOverlapZone
            nOverlap = nOverlap + 1;
            vMatch = t==t(n);
            vMatch(n) = false;
            inds = find(vMatch);
            if isempty(inds)
                matchIdx = NaN;
            else
                matchIdx = inds(1);
            end
            iOverlap(nOverlap, 1) = n;
            iMatch(nOverlap, 1) = matchIdx;
        else
            lastValidTime = t(n);
        end
    end
end

function data = fixOverlappingFrames(data, iOverlap, iMatch, mode)
    % New approach: treat overlapping frames as either duplicated or 
    % out-of-order. Possibilities are (1) discard all frames from time
    % regions were identified as ovelapping, or (2) discard only frames
    % which have exactly duplicated timestamps, and sort all remaining
    % frames in timestamp order.
    t = data.frameTimestamp;
    
    if strcmpi(mode, 'sort')
        [~, iSort] = sort(t);
        vKeep = true(size(t));
        vKeep(iOverlap(~isnan(iMatch))) = false;
        iSort = iSort(vKeep(iSort));
        frameInds = iSort;
    elseif strcmpi(mode, 'discard')
        frameInds = (1:numel(t));
        frameInds(iOverlap) = [];
    end
    
    fds = fieldnames(data);
    for f = 1:numel(fds)
        fd = fds{f};
        data.(fd) = data.(fd)(frameInds, :);
    end
    t = data.frameTimestamp;
    assert(all(diff(t)>0));
end

function data = fixTimestamps(data, P)
% Ensures timestamps are monotonically ascending

t = data.frameTimestamp;

increment = P.fixTimestampsIncrement;
thresh = P.fixTimestampsThresholdFrames;

dt = diff(t);
iNeg = find(dt<0);
iNegBef = max(1, iNeg-thresh+1);
iNegAft = min(numel(t), iNeg+thresh);

iBad = [];
for n = 1:numel(iNeg)
    indsBef = iNegBef(n):iNeg(n);
    indsAft = (iNeg(n)+1 : iNegAft(n));
    if min(t0(indsBef)) > max(t(indsAft))
        iBad(end+1, 1) = iNeg(n);
    end
end

nBad = numel(iBad);
nFrames = numel(t);

if nBad > 0
%     fprintf('Timestamps are not monotonically ascending. Fixing %u negative jumps.\n', nBad);
    for i = 1:nBad
        inds = (iBad(i)+1 : nFrames)';
        tShift = t(inds(1)-1) - t(inds(1)) + increment;
        t(inds) = t(inds) + tShift;
    end
end

data.frameTimestamp = t;

end