function splitMtvFile(sourceFile, targetFiles)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

import matmot.FormatSpec

CHUNK_SIZE = 2^15;

[pth, fn, ext] = fileparts(sourceFile);
rootDir = pth;
metaDir = fullfile(pth, [fn '_meta']);
if ~exist(metaDir, 'dir')
    error('Pre-merged .meta directory "%s" not found.', metaDir);
end
tmp = dir(fullfile(metaDir, '*.meta'));
metaFiles = {tmp.name};
nFiles = numel(metaFiles);
for f = 1:nFiles
    fpth = fullfile(metaDir, metaFiles{f});
    [pth, fn, ext] = fileparts(fpth);
    txt = FormatSpec.readMetaText(fpth);
    meta = FormatSpec.parseMetaText(txt);
    fileSizes(f) = meta.file_size;
    fileInds(f) = str2num(fn);
end

[fileInds, iSort] = sort(fileInds);
assert(isequal(fileInds, 1:nFiles));
metaFiles = metaFiles(iSort);
fileSizes = fileSizes(iSort);
splitAtBytes = cumsum(fileSizes);

nSubfiles = numel(splitAtBytes);

if nargin < 2 || isempty(targetFiles)
    targetFiles = arrayfun(@(idx) {sprintf('%u.meta', idx)}, fileInds);
end

txt = FormatSpec.readMetaText(sourceFile);
srcMeta = FormatSpec.parseMetaText(txt);
fidIn = fopen(sourceFile, 'r');

for f = 1:nSubfiles
    srcMetaFile = fullfile(metaDir, metaFiles{f});
    
    % Write output files to same dir as source unless specified
    [pth, fn, ext] = fileparts(targetFiles{f});
    if isempty(pth)
        pth = rootDir;
    end
    targetFile = fullfile(pth, [fn '.mtv']);
    targetMetaFile = fullfile(pth, [fn '.meta']);
    txt = FormatSpec.readMetaText(srcMetaFile);
    meta = FormatSpec.parseMetaText(txt);
    FormatSpec.writeMeta(targetMetaFile, FormatSpec.VERSION, meta);
    fidOut = fopen(targetFile, 'w');
    nextSplit = splitAtBytes(f);
    
    while ~feof(fidIn) && ftell(fidIn) < nextSplit
        buff = fread(fidIn, CHUNK_SIZE, '*uint8');
        splitOffset = ftell(fidIn) - nextSplit;
        if splitOffset > 0
            buff = buff(1: (CHUNK_SIZE-splitOffset));
            fseek(fidIn, nextSplit, 'bof');
        end
        fwrite(fidOut, buff);
    end
    fclose(fidOut);
    
end

fclose(fidIn);