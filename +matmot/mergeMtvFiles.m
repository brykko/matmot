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
inp.parse(varargin{:});
P = inp.Results;

CHUNK_N_FRAMES = 1024;

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
[pth, fn, ext] = fileparts(mergedFilename);
mergedFilenameMeta = fullfile(pth, [fn '.meta']);
mergedMetaDir = fullfile(pth, [fn '_meta']);
if exist(mergedMetaDir, 'dir')
    rmdir(mergedMetaDir, 's');
end
mkdir(mergedMetaDir);

nBytesTotal = 0;
fidIn = 0;

try
    nFramesTotal = 0;
    
    for f = 1:nFiles
        fn = filenames{f};
        [meta, fileVer, offset] = readMeta(fn);
        nBytesFrame = bytesPerFrame(meta.n_rigid_bodies, meta.n_markers);
        fidIn = fopen(fn, 'r');
        fseek(fidIn, offset, 'bof');
        
        % Copy data to new file
        fprintf('Writing data from source file "%s" to target file "%s"...', fn, mergedFilename);
        nBytes = 0;
        nFiles = 0;
        
        while nFiles < meta.n_frames
            bytes = fread(fidIn, [nBytesFrame, CHUNK_N_FRAMES], '*uint8');
            nFiles = nFiles + size(bytes, 2);
            if nFiles > meta.n_frames
                nDrop = nFiles - meta.n_frames;
                bytes = bytes(:, 1:end-nDrop);
            end
            fwrite(fidOut, bytes);
            nBytes = nBytes + numel(bytes);
        end
        
        nBytesTotal = nBytesTotal + nBytes;
        fclose(fidIn);
        fprintf('done.\n');
        
        % Write meta text to subdir
        meta.file_size = nBytes;
        fnMeta = fullfile(mergedMetaDir, sprintf('%u.meta', f));
        writeMeta(fnMeta, fileVer, meta)
        nFramesTotal = nFramesTotal + meta.n_frames;
    end
    
    % Write .meta for target file
    meta.n_frames = nFramesTotal;
    meta.file_size = nBytesTotal;
    meta.composite_file = 'true';
    writeMeta(mergedFilenameMeta, Version(), meta)
    
    fclose(fidOut);
    fprintf('All files written.\n');
    succeeded = true;
    msg = '';
    
catch e
    % If anything fails, close everything and about, reporting failure
    warning('Merging files failed: "%s". Aborting.', e.message);
    fids = [fidIn, fidOut];
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