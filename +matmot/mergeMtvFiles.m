function mergeMtvFiles(mergedFilename, filenames)
%MERGEMTVFILES concatenate multiple .mtv files
%
% This function concatenates the data from a series of .mtv files. All the
% files must contain the same number of rigid bodies and markers.
%
% MERGEMTVFILES(FNSAVE, FN1, FN2, FN3 ...) concatenates the data from .mtv
% files FN1, FN2, FN3 (etc.), saving the joined data in file FNSAVE. FNSAVE
% will always take the file header from the first input filename.

CHUNK_SIZE = 2^15;

% Check all files exist
% filenames = varargin;
nFiles = numel(filenames);
for f = 1:nFiles
    fn = filenames{f};
    assert(exist(fn, 'file') > 0, 'File "%s" could not be found', fn);
end

fidOut = fopen(mergedFilename, 'W');
headerLen = matmot.FormatSpec.HEADER_LENGTH;

fidIn = fopen(filenames{1}, 'r');
headerBytes = fread(fidIn, headerLen, '*uint8');
fwrite(fidOut, headerBytes);
fclose(fidIn);

for f = 1:nFiles
    fn = filenames{f};
    fprintf('Writing file "%s"...', fn);
    fidIn = fopen(fn, 'r');
    fseek(fidIn, headerLen, 'bof');
    while ~feof(fidIn)
        bytes = fread(fidIn, CHUNK_SIZE, '*uint8');
        fwrite(fidOut, bytes);
    end
    fclose(fidIn);
    fprintf('done\n', fn);
end

fclose(fidOut);
fprintf('All files written.\n');

end

