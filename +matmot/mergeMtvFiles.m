function [succeeded, msg] = mergeMtvFiles(mergedFilename, filenames, varargin)
%MERGEMTVFILES concatenate multiple .mtv files
%
% This function concatenates the data from a series of .mtv files. All the
% files must contain the same number of rigid bodies and markers.
%
% MERGEMTVFILES(FNSAVE, FN1, FN2, FN3 ...) concatenates the data from .mtv
% files FN1, FN2, FN3 (etc.), saving the joined data in file FNSAVE. FNSAVE
% will always take the file header from the first input filename.

import matmot.FormatSpec.*

inp = inputParser();
inp.addParameter('force', false);
inp.addParameter('verbose', false);
inp.addParameter('resetTimestamps', true);
inp.parse(varargin{:});
P = inp.Results;

CHUNK_N_FRAMES = 1024;
FILE_SPACING_TIME = 1;
FILE_SPACING_IDX = 1;

% Check all files exist
nFiles = numel(filenames);
for f = 1:nFiles
    fn = filenames{f};
    assert(exist(fn, 'file') > 0, 'File "%s" could not be found', fn);
end

if exist(mergedFilename, 'file') && ~P.force
    warning('File "%s" already exists, aborting.', mergedFilename);
    return;
end

fidOut = fopen(mergedFilename, 'W');
[pth, fn] = fileparts(mergedFilename);
mergedFilenameMeta = fullfile(pth, [fn '.meta']);
mergedMetaDir = fullfile(pth, [fn '_meta']);

if exist(mergedMetaDir, 'dir')
    rmdir(mergedMetaDir, 's');
end
mkdir(mergedMetaDir);

% Initialize output info text file
mergedInfoFilename = fullfile(pth, [fn '.mergeinfo']);
fidInf = fopen(mergedInfoFilename, 'w');
fprintf(fidInf, 'reset_timestamps=%u\r\n', P.resetTimestamps);
fprintf(fidInf, 'output_file="%s"\r\n', mergedFilename);
fprintf(fidInf, 'file_spacing_time=%.6f\r\n', FILE_SPACING_TIME);
fprintf(fidInf, 'file_spacing_index=%u\r\n', FILE_SPACING_IDX);

nBytesTotal = 0;
fidIn = 0;
nFramesTotal = 0;
lastTimestamp = 0;
lastIdx = int32(0);

try
    for f = 1:nFiles
        
        fn = filenames{f};
        [meta, fileVer, offset] = readMeta(fn);
        nBytesFrame = bytesPerFrame(meta.n_rigid_bodies, meta.n_markers);
        
        fidIn = fopen(fn, 'r');
        fseek(fidIn, offset, 'bof');
        
        % Copy data to new file
        fprintf('Writing data from source file "%s" to target file "%s"...\n', fn, mergedFilename);
        nBytes = 0;
        nFrames = 0;
        nBadTimes = 0;
        nBadInds = 0;
        firstChunk = true;
        fileBadFrames = false;
        
        while nFrames < meta.n_frames
            bytes = fread(fidIn, [nBytesFrame, CHUNK_N_FRAMES], '*uint8');
            bytesFrameI = bytes(1:4, :);
            bytesFrameT = bytes((5:12), :);
            tFrame = typecast(bytesFrameT(:), 'double');
            iFrame = typecast(bytesFrameI(:), 'int32');
            if firstChunk && P.resetTimestamps
                gapT = lastTimestamp - tFrame(1) + FILE_SPACING_TIME;
                gapI = lastIdx - iFrame(1) + FILE_SPACING_IDX;
            end
            tFrame = tFrame+gapT;
            iFrame = iFrame+gapI;
            
            nBadInds = nBadInds + numel(find(iFrame < lastIdx));
            nBadTimes = nBadTimes + numel(find(tFrame < lastTimestamp));
            
            bytes(1:4, :) = reshape(typecast(iFrame, 'uint8'), 4, []);
            bytes((5:12), :) = reshape(typecast(tFrame, 'uint8'), 8, []);
            
            nFrames = nFrames + size(bytes, 2);
            if nFrames > meta.n_frames
                nDrop = nFrames - meta.n_frames;
                bytes = bytes(:, 1:end-nDrop);
            end
            fwrite(fidOut, bytes);
            nBytes = nBytes + numel(bytes);
            
            firstChunk = false;
            lastTimestamp = tFrame(end);
            lastIdx = iFrame(end);
            
        end
        
        fileFrameRange = [nFramesTotal+1, nFramesTotal + meta.n_frames];
        fprintf( ...
            fidInf, ...
            'file_index=%u, name="%s", composite_frame_inds=(%u-%u), t_adjust=%.6f, i_adjust=%d, n_bad_frames=%u, n_bad_inds=%u\r\n', ...
            f, fn, fileFrameRange(1), fileFrameRange(2), gapT, gapI, nBadTimes, nBadInds);
        
        nBytesTotal = nBytesTotal + nBytes;
        fclose(fidIn);
        fprintf('done.\n');
        
        if nBadTimes || nBadInds
            warning('matmot:mergeMtvFiles:badSourceFrames', ...
                'Input data file %s appears to contain bad (possibly duplicated) frames: %u timestamps, %u indices were out of order.', fn, ...
                nBadTimes, nBadInds);
        end
        
        % Write meta text to subdir
        meta.file_size = nBytes;
        fnMeta = fullfile(mergedMetaDir, sprintf('%u.meta', f));
        writeMeta(fnMeta, fileVer, meta)
        nFramesTotal = fileFrameRange(2);
    end
    
    % Write .meta for target file
    meta.n_frames = nFramesTotal;
    meta.file_size = nBytesTotal;
    meta.composite_file = 'true';
    meta.composite_version = VERSION;
    writeMeta(mergedFilenameMeta, Version(), meta)
    
    fclose(fidOut);
    fclose(fidInf);
    fprintf('All files written.\n');
    succeeded = true;
    msg = '';
    
catch e
    % If anything fails, close everything and about, reporting failure
    warning('Merging files failed: "%s". Aborting.', e.message);
    fids = [fidIn, fidOut, fidInf];
    for f = 1:2
        try
            fclose(fids(f));
        catch
        end
    end
    succeeded = false;
    msg = e.message;
end

end