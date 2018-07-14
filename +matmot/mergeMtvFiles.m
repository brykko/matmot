function [succeeded, msg] = mergeMtvFiles(mergedFilename, filenames, varargin)
%MERGEMTVFILES concatenate multiple .mtv files
%
% This function concatenates the data from a series of .mtv files. All the
% files must contain the same number of rigid bodies and markers.
%
% MERGEMTVFILES(FNSAVE, FN1, FN2, FN3 ...) concatenates the data from .mtv
% files FN1, FN2, FN3 (etc.), saving the joined data in file FNSAVE. FNSAVE
% will always take the file header from the first input filename.

import matmot.FormatSpec

inp = inputParser();
inp.addParameter('force', false);
inp.parse(varargin{:});
P = inp.Results;

CHUNK_SIZE = 2^15;

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
mergedMetaDir = [fn '_meta'];
if exist(mergedMetaDir, 'dir')
    rmdir(mergedMetaDir, 's');
end
mkdir(mergedMetaDir);

nBytesTotal = 0;

try
    nFramesTotal = 0;
    
    for f = 1:nFiles
        fn = filenames{f};
        metaTxt = FormatSpec.readMetaText(fn);
        meta = FormatSpec.parseMetaText(metaTxt);
        fileVer = meta.matmot_version;
        useMetaFile = fileVer(2) > 1;
        
        fprintf('Writing file "%s"...\n', fn);
        fidIn = fopen(fn, 'r');
        if ~useMetaFile
            fseek(fidIn, FormatSpec.HEADER_LENGTH, 'bof');
        end
        
        % Copy data to new file
        nBytes = 0;
        while ~feof(fidIn)
            bytes = fread(fidIn, CHUNK_SIZE, '*uint8');
            fwrite(fidOut, bytes);
            nBytes = nBytes + numel(bytes);
        end
        nBytesTotal = nBytesTotal + nBytes;
        fclose(fidIn);
        
        % Write meta text to subdir
        if fileVer(2) < 2
            meta.n_frames = FormatSpec.fileNFrames(fn);
            meta.file_size = nBytes;
        end
        
        fnMeta = fullfile(mergedMetaDir, sprintf('%u.meta', f));
        FormatSpec.writeMeta(fnMeta, fileVer, meta)
        nFramesTotal = nFramesTotal + meta.n_frames;
    end
    
    % Write .meta for target file
    meta.n_frames = nFramesTotal;
    meta.file_size = nBytesTotal;
    meta.composite_file = 'true';
    FormatSpec.writeMeta(mergedFilenameMeta, FormatSpec.VERSION, meta)
    
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