function [inds, times, data] = extractSyncPulses(filepath, varargin)
%EXTRACTSYNCPULSES extracts IR-LED sync pulse times from Motive data
%
%   This function can be used to detect the locations and pulse times of a
%   set of IR-LEDs positioned in the recorded volume. All LEDs must be
%   pulsed on and off in unison, since the detection algorithm uses the
%   correlations of voxel occupancies to detect which voxels the LEDs are
%   located in.
%
%   [INDS, TIMES] = EXTRACTSYNCPULSES(FILEPATH) processes the file
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

inp = inputParser();
inp.addParameter('plot', false);
inp.addParameter('figure', []);
inp.addParameter('nCheck', 10);
inp.addParameter('voxelCorrThreshold', 0.9);
inp.addParameter('nLedsMin', 2);
inp.addParameter('nLedsMax', 4);
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

% Calculate voxel/pixel-occupancy
mx = data.mx(:);
my = data.my(:);
mz = data.mz(:);
mPos = {data.mx, data.my, data.mz};
nSamples = size(data.mx, 1);

[histCounts, histEdges, histCenters] = matmot.external.histcn([mx, my, mz]);
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
    [idx(1), idx(2), idx(3)] = ind2sub(size(histCounts), indsSort(b));
    hiOccupancyVoxelInds(b, :) = idx;
    inBin = true(nSamples, 1);
    for dim = 1:3
        inBin = inBin & anyMarkerInBin(histEdges{dim}, idx(dim), mPos{dim});
    end
    if P.plot
        x = data.frameTimestamp;
        y = inBin/2 + b;
        line(ax, x, y, 'color', 'k');
    end
    voxelOccupancies(:, b) = inBin;
end

r = corr(voxelOccupancies);
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
    combs = nchoosek(1:P.nCheck, n);
    nCombs = size(combs, 1);
    for c = 1:nCombs
        % Get the subset of voxels
        comb = combs(c, :);
        % Calculate all pairwise combinations in the voxel subset
        combsr = nchoosek(1:n, 2);
        i = comb(combsr(:, 1));
        j = comb(combsr(:, 2));
        idx = sub2ind(size(r), i, j);
        % Get the minimum correlation between any pair
        minCorr = min(r(idx));
    
        % As soon as all pairs are correlated above the threshold, assume
        % these voxels all contain the LEDs.
        if minCorr > P.voxelCorrThreshold
            foundLedVoxels = true;
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
    nLedVoxels, minCorr);
fprintf('LED voxel coordinates (x, y, z):\n');
for n = 1:nLedVoxels
    fprintf('(%.3f, %.3f, %.3f), ', ...
        ledVoxelCoords(n, 1), ...
        ledVoxelCoords(n, 2), ...
        ledVoxelCoords(n, 3));
end
fprintf('\n\n');
if P.plot
   ax = subplot(sply, splx, 5, 'color', 'k');
   pts = ledVoxelCoords([1:end, 1], :);
   line(ax, pts(:, 1), pts(:, 2), pts(:, 3), 'marker', '.', 'markerSize', 20, ...
       'color', 'r', 'lineStyle', 'none');
   pts = data.pos;
   % Plot the rigid body path. Need to reverse x-direction to make it match
   % marker coords (is this correct?)
   line(ax, -pts(:, 1), pts(:, 2), pts(:, 3), 'color', 'w', 'lineWidth', 0.5);
   axis(ax, 'equal');
   xlabel(ax, 'x');
   ylabel(ax, 'y');
   zlabel(ax, 'z');
   title(ax, 'Identified LED voxel locations');
end

pulseStates = any(voxelOccupancies(:, ledVoxelInds), 2);

% Extract the pulse onset and offset times
p1 = pulseStates(2:end);
p0 = pulseStates(1:end-1);
isU = [false; p1 & ~p0];
isD = [false; ~p1 & p0];

iU = find(isU);
iD = find(isD);

if iU(1) > iD(1)
    iD(1) = [];
end

if iU(end) > iD(end)
    iU(end) = [];
end

inds = [iU iD];
times = data.frameTimestamp(inds);

pulseLen = diff(times, [], 2);
pulseInterval = times(2:end, 1) - times(1:end-1, 2);
fprintf('Median pulse length = %.1f ms, median interpulse interval = %.1f ms\n', ...
    median(pulseLen)*1e3, median(pulseInterval)*1e3);

end

function v = anyMarkerInBin(binEdges, binIdx, positions)
    bin = binEdges(binIdx + [0, 1]);
    v = any(positions > bin(1) & positions <= bin(2), 2);
end