function nbytes = fileSize(filename)
fid = fopen(filename, 'r');
fseek(fid, 0, 'eof');
nbytes = ftell(fid);
fclose(fid);
end

