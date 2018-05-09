function [events, data] = extractSyncPulses(filepath, varargin)
%EXTRACTSYNCPULSES extracts IR-LED sync pulse times from Motive data
%
%   This function can be used to detect the locations and pulse times of a
%   set of IR-LEDs positioned in the recorded volume. All LEDs must be
%   pulsed on and off in unison, since the detection algorithm uses the
%   correlations of voxel occupancies to detect which voxels the LEDs are
%   located in.
%
%   [EVENTS, TIMES] = EXTRACTSYNCPULSES(FILEPATH) processes the file
%   containing streamed Motive data specified by FILEPATH. The returned
%   variables INDS and TIMES are two-column matrices indicating the pulse 
%   onset and offset frame indices and times, respectively.
%
%   [...] = EXTRACTSYNCPULSES(FILEPATH, PRM, VAL, ...) specifies optional 
%   parameter/value pairs. Available parameters are as follows:
%   
%   'plot' (default FALSE): generate a summary figure.
%
%   'figure' (default gcf): if parameter 'plot' is enabled, this 
%       specifies which figure will be used for plotting.
%
%   'nCheck' (default 10): the number of highest-occupancy voxels that 
%       will be checked.
%
%   'voxelCorrThreshold' (default 0.9): defines the minimum correlation
%       between a subset of the highest-occupancy voxels that qualifies 
%       them as containing LEDs.
%
%   'nLedsMin'  (default 2) the minimum number of LEDs to look for. This 
%       can be set equal to the number of LEDs present in the
%       recording, or fewer to allow for the possibility that some LEDs may 
%       not have been detected.
%
%   'nLedsMax' (integer, default 4) the maximum number of LEDs to look for.
%       This should be equal to the number of LEDs used in the recording.
%
%   'validPositionRange' (3-element cell array, default 
%       {[-5, 5], [-5, 5], [-5, 5]} ) x, y, z coordinates of the limits for
%       valid positions. If any marker positions are outside this range,
%       they will be set to NaN and will not be registered.

inp = inputParser();
inp.addParameter('plot', false);
inp.addParameter('figure', []);
inp.addParameter('nCheck', 50);
inp.addParameter('voxelCorrThreshold', 0.9);
inp.addParameter('nLedsMin', 3);
inp.addParameter('nLedsMax', 4);
inp.addParameter('ledMinPercentOccupancy', 20);
inp.addParameter('validPositionRange', {[-5 5], [-5 5], [-5 5]});
inp.parse(varargin{:});
P = inp.Results;

if P.plot
    if isempty(P.figure)
        fig = gcf();
    else
        fig = P.figure;
    end
    splx = 3;
    sply = 2;
end

data = matmot.loadMtvFile(filepath);
nFrames = size(data.rbx, 1);

% Concatenate all marker positions and generate 3D histogram to determine
% where the highest-occupancy voxels are. Discard any spurious positions
% that lie outside the valid range.
mPos = {data.mx, data.my, data.mz};
for n = 1:3
    tmp = mPos{n};
    rng = P.validPositionRange{n};
    validPos = tmp >= rng(1) & tmp <= rng(2);
    tmp(~validPos) = nan;
    mPos{n} = tmp;
end
    
nSamples = size(data.mx, 1);
[histCounts, histEdges, histCenters] = matmot.external.histcn( ...
    [mPos{1}(:), mPos{2}(:), mPos{3}(:)]);
histCounts2d = squeeze(sum(histCounts, 2));

if P.plot
    % voxel occupancy histogram
    ax = subplot(sply, splx, 1, 'parent', fig, 'color', 'k');
    alpha = histCounts / max(histCounts(:));
    cdata = ones([size(histCounts) 3]);
    matmot.external.vol3d('cdata', cdata, 'alpha', alpha, 'parent', ax);
    xlabel(ax, 'x');
    ylabel(ax, 'y');
    zlabel(ax, 'z');
    axis(ax, 'equal');
    title('Marker distribution (3d)')
    ax.CLim = [0 1];
    colormap(ax, 'gray');
    colorbar(ax);
    axis(ax, 'image');
    
    % x/z pixel occupancy histogram
    ax = subplot(sply, splx, 2, 'parent', fig);
    colormap(ax, gray());
    cdata = histCounts2d / max(histCounts2d(:));
    imagesc(ax, histCenters{1}, histCenters{3}, cdata);
    colorbar(ax);
    title(ax, 'Marker distribution (horz. plane)')
    axis(ax, 'xy');
    xlabel(ax, 'x');
    ylabel(ax, 'z');
    axis(ax, 'image');
end

% For each high-density bin, find its occupancy onset and offset times
[~, indsSort] = sort(histCounts(:), 'descend');

if P.plot
    ax = subplot(sply, splx, 3, 'parent', fig);
    ax.YTick = 1:P.nCheck;
    title(ax, 'Highest-density voxel occupancy traces');
    ylabel(ax, 'Voxel #');
    xlabel(ax, 'Time');
    ax.XLim = data.frameTimestamp([1, end]);
end

voxelOccupancies = false(nSamples, P.nCheck);

for b = 1:P.nCheck
    % Get the edges of the current voxel
    [inds(1), inds(2), inds(3)] = ind2sub(size(histCounts), indsSort(b));
    hiOccupancyVoxelInds(b, :) = inds;
    inBin = true(nSamples, 1);
    for dim = 1:3
        inBin = inBin & anyMarkerInBin(histEdges{dim}, inds(dim), mPos{dim});
    end
    if P.plot
        x = data.frameTimestamp;
        y = inBin/2 + b;
        line(ax, x, y, 'color', 'k');
    end
    voxelOccupancies(:, b) = inBin;
end

voxelPercentOccupancies = 100 * sum(voxelOccupancies) / nFrames;
validVoxels = voxelPercentOccupancies >= P.ledMinPercentOccupancy;
hiOccupancyVoxelInds = hiOccupancyVoxelInds(validVoxels, :);
voxelOccupancies = voxelOccupancies(:, validVoxels);
r = corr(voxelOccupancies);
nCheck = numel(find(validVoxels));

if P.plot
    ax = subplot(sply, splx, 4, 'parent', fig);
    imagesc(ax, r);
    ax.XTick = 1:P.nCheck;
    ax.YTick = 1:P.nCheck;
    colormap(ax, hot());
    cb = colorbar(ax);
    ylabel(cb, 'r');
    axis(ax, 'image');
    title(ax, 'Voxel occupancy correlation')
end

% Find the voxels containing the LEDs
foundLedVoxels = false;

for n = P.nLedsMax : -1 : P.nLedsMin
    combs = nchoosek(1:nCheck, n);
    nCombs = size(combs, 1);
%     maxCorr = zeros(nCombs, 1);
    minCorr = zeros(nCombs, 1);
    for c = 1:nCombs
        % Get the subset of voxels
        comb = combs(c, :);
        % Calculate all pairwise combinations in the voxel subset
        combsr = nchoosek(1:n, 2);
        i = comb(combsr(:, 1));
        j = comb(combsr(:, 2));
        inds = sub2ind(size(r), i, j);
        % Get the minimum correlation between any pair
%         maxCorr(c) = max(r(inds));
        minCorr(c) = min(r(inds));
    end
    
        % As soon as all pairs are correlated above the threshold, assume
        % these voxels all contain the LEDs.
        if any(minCorr > P.voxelCorrThreshold)
            foundLedVoxels = true;
            [~, idx] = max(minCorr);
            comb = combs(idx, :);
            bestMinCorr = minCorr(idx);
            ledVoxelInds = comb;
            nLedVoxels = numel(comb);
            for nn = 1:n
                ledVoxelBins(nn, :) = hiOccupancyVoxelInds(comb(nn), :);
                for dim = 1:3
                    ledVoxelCoords(nn, dim) = histCenters{dim}(ledVoxelBins(nn, dim));
                end
            end
            break;
        end
    
    if foundLedVoxels
        break;
    end
    
end

if ~foundLedVoxels
    error('Failed to find group of %u-%u voxels with r > %.2f', ...
        P.nLedsMin, P.nLedsMax, P.voxelCorrThreshold)
end

fprintf('Found %u voxels with minimum r = %.2f\n\n', ...
    nLedVoxels, bestMinCorr);
fprintf('LED voxel coordinates (x, y, z):\n');
for n = 1:nLedVoxels
    fprintf('(%.3f, %.3f, %.3f), ', ...
        ledVoxelCoords(n, 1), ...
        ledVoxelCoords(n, 2), ...
        ledVoxelCoords(n, 3));
end
fprintf('\n\n');

if P.plot
   clear h
   ax = subplot(sply, splx, 5, 'color', 'k');
   pts = ledVoxelCoords([1:end, 1], :);
   h(1) = line(ax, pts(:, 1), pts(:, 2), pts(:, 3), 'marker', 'o', 'markerSize', 5, ...
       'color', 'r', 'lineStyle', 'none');
   % Plot the rigid body path
   h(2) = line(ax, data.rbx, data.rby, data.rbz, 'color', 'w', 'lineWidth', 0.5);
   axis(ax, 'equal');
   xlabel(ax, 'x');
   ylabel(ax, 'y');
   zlabel(ax, 'z');
   title(ax, 'Identified LED voxel locations');
   axis(ax, 'equal');
   rotate3d(ax);
end

pulseStates = any(voxelOccupancies(:, ledVoxelInds), 2);

% Extract the pulse onset and offset times
events = detectEvents(pulseStates, data.frameTimestamp);
% p1 = pulseStates(2:end);
% p0 = pulseStates(1:end-1);
% isU = [false; p1 & ~p0];
% isD = [false; ~p1 & p0];
% 
% iU = find(isU);
% iD = find(isD);
% 
% if iU(1) > iD(1)
%     iD(1) = [];
% end
% 
% if iU(end) > iD(end)
%     iU(end) = [];
% end

% inds = [[events.start]' [events.stop]'];

% Convert frame indices to times (only for verbose feedback)
% times = data.frameTimestamp(inds);
% pulseLen = diff(times, [], 2);
pulseLen = [events.tStop]-[events.tStart];
pulseInterval = [events(2:end).tStart] - [events(1:end-1).tStop];
% pulseInterval = times(2:end, 1) - times(1:end-1, 2);
fprintf('Median pulse length = %.1f ms, median interpulse interval = %.1f ms\n', ...
    median(pulseLen)*1e3, median(pulseInterval)*1e3);

end

function v = anyMarkerInBin(binEdges, binIdx, positions)
    nBins = numel(binEdges);
    if binIdx < nBins
        bin = binEdges(binIdx + [0, 1]);
        v = positions > bin(1) & positions <= bin(2);
    else
        v = positions == binEdges(end);
    end
    v = any(v, 2);
end