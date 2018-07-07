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
% filenames = varargin;
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
filepathSum = [mergedFilename '.cat'];
fidSum = fopen(filepathSum, 'w');
fidIn = fopen(filenames{1}, 'r');

try
    
    headerBytes = fread(fidIn, FormatSpec.HEADER_LENGTH, '*uint8');
    fwrite(fidOut, headerBytes);
    fclose(fidIn);
    
    for f = 1:nFiles
        fn = filenames{f};
        fprintf('Writing file "%s"...\n', fn);
        fidIn = fopen(fn, 'r');
        fseek(fidIn, FormatSpec.HEADER_LENGTH, 'bof');
        nBytesData = 0;
        while ~feof(fidIn)
            bytes = fread(fidIn, CHUNK_SIZE, '*uint8');
            fwrite(fidOut, bytes);
            nBytesData = nBytesData + numel(bytes);
        end
        fclose(fidIn);
        fprintf(fidSum, 'file_index=%u, n_bytes=%u, path="%s"\r\n', ...
            f, nBytesData, filenames{f});
    end
    
    fclose(fidOut);
    fclose(fidSum);
    fprintf('All files written.\n');
    succeeded = true;
    msg = '';
    
catch e
    % If anything fails, close everything and about, reporting failure
    warning('Merging files failed: "%s". Aborting.', e.message);
    fids = [fidSum, fidIn, fidOut];
    for f = 1:3
        try
            fclose(fids(f));
        catch
        end
    end
    succeeded = false;
    msg = e.message;
end

end