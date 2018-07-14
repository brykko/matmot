classdef FormatSpec < handle
    %FormatSpec file format-related constants and functions
    
    properties (Constant)
        
        FRAME_FIELD_NAMES = {
            'frameIdx'
            'frameTimestamp'
            'frameLatency'}
        
        FRAME_ENCODINGS = {
            'int32'
            'double'
            'single'}
        
        RB_FIELD_NAMES = {
            'rbx'
            'rby'
            'rbz'
            'rbqx'
            'rbqy'
            'rbqz'
            'rbqw'
            'rbError'
            'rbTracked'}
        
        RB_ENCODINGS = {
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'uint8'}
        
        MARKER_FIELD_NAMES = {
            'mx'
            'my'
            'mz'
            'msz'
            'mres'}
        
        MARKER_ENCODINGS = {
            'single'
            'single'
            'single'
            'single'
            'single'}
        
        HEADER_LENGTH = 2^14
        
        TIMESTAMP_JOIN_INCR = 10000;
        
classdef FormatSpec < handle
    %FormatSpec file format-related constants and functions
    
    properties (Constant)
        
        FRAME_FIELD_NAMES = {
            'frameIdx'
            'frameTimestamp'
            'frameLatency'}
        
        FRAME_ENCODINGS = {
            'int32'
            'double'
            'single'}
        
        RB_FIELD_NAMES = {
            'rbx'
            'rby'
            'rbz'
            'rbqx'
            'rbqy'
            'rbqz'
            'rbqw'
            'rbError'
            'rbTracked'}
        
        RB_ENCODINGS = {
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'uint8'}
        
        MARKER_FIELD_NAMES = {
            'mx'
            'my'
            'mz'
            'msz'
            'mres'}
        
        MARKER_ENCODINGS = {
            'single'
            'single'
            'single'
            'single'
            'single'}
        
        HEADER_LENGTH = 2^14
        
        VERSION = '0.2.0';
        
    end
    
    methods (Static)
        
        function S = rbFields()
            % Field encodings as a struct array
            import matmot.FormatSpec
            [byteInds, nBytes] = FormatSpec.getByteInds(FormatSpec.RB_ENCODINGS);
            S = struct( ...
                'name', FormatSpec.RB_FIELD_NAMES, ...
                'byte_inds', num2cell(byteInds'), ...
                'n_bytes', num2cell(nBytes'), ...
                'encoding', FormatSpec.RB_ENCODINGS);
        end
        
        function S = markerFields()
            % Marker field encodings as a struct array
            import matmot.FormatSpec
            [byteInds, nBytes] = FormatSpec.getByteInds(FormatSpec.MARKER_ENCODINGS);
            S = struct( ...
                'name', FormatSpec.MARKER_FIELD_NAMES, ...
                'byte_inds', num2cell(byteInds'), ...
                'n_bytes', num2cell(nBytes'), ...
                'encoding', FormatSpec.MARKER_ENCODINGS);
        end
        
        function S = basicFields()
            % Marker field encodings as a struct array
            import matmot.FormatSpec
            [byteInds, nBytes] = FormatSpec.getByteInds(FormatSpec.FRAME_ENCODINGS);
            S = struct( ...
                'name', FormatSpec.FRAME_FIELD_NAMES, ...
                'byte_inds', num2cell(byteInds'), ...
                'n_bytes', num2cell(nBytes'), ...
                'encoding', FormatSpec.FRAME_ENCODINGS);
        end
        
        function nBytes = bytesBasic()
            % Number of bytes for basic frame fields
            import matmot.FormatSpec
            fields = FormatSpec.basicFields();
            nBytes = sum([fields.n_bytes]);
        end
        
        function nBytes = bytesPerRb()
            % Number of bytes needed to encode data of rigid body
            import matmot.FormatSpec
            fields = FormatSpec.rbFields();
            nBytes = sum([fields.n_bytes]);
        end
        
        function nBytes = bytesPerMarker()
            % Number of bytes needed to encode data of one marker
            import matmot.FormatSpec
            fields = FormatSpec.markerFields();
            nBytes = sum([fields.n_bytes]);
        end
        
        function [encoding, n] = convertEncoding(encoding)
            % Decompose a composite encoding string e.g. '4*single'
            if contains(encoding, '*')
                idx = find(encoding=='*');
                n = str2double(encoding(1:idx-1));
                encoding = encoding(idx+1:end);
            else
                n = 1;
            end
        end
        
        function nBytes = encodingNBytes(encoding)
            % Number of bytes corresponding to a composite encoding
            [encoding, n] = matmot.FormatSpec.convertEncoding(encoding);
            x = zeros(1, n, encoding);
            nBytes = numel(typecast(x, 'uint8'));
        end
        
        function [inds, nBytes] = getByteInds(encodings)
            % Byte indices and numbers
            import matmot.FormatSpec
            idx0 = 1;
            for n = 1:numel(encodings)
                encoding = encodings{n};
                nBytes(n) = FormatSpec.encodingNBytes(encoding);
                inds(n) = idx0;
                idx0 = idx0 + nBytes(n);
            end
        end
        
        function nBytes = bytesPerFrame(nRbs, nMarkers)
            %BYTESPERFRAME number of bytes needed to encode one frame.
            %
            % NBYTES = BYTESPERFRAME(NRBS, NMARKERS) returns the number of
            % bytes NBYTES needed to encode a frame containing the rigid
            % body count NRBS and the marker count NMARKERS.
            import matmot.FormatSpec
            nBytes = ...
                FormatSpec.bytesBasic() + ...
                nRbs*FormatSpec.bytesPerRb() + ...
                nMarkers*FormatSpec.bytesPerMarker();
        end
        
        function writeMeta(filename, version, params)
            paramNames = fieldnames(params);
            fid = fopen(filename, 'w');
            if ~ischar(version)
                version = sprintf('%u.%u.%u', version(1), version(2), version(3));
            end
            fprintf(fid, 'MatMot Streamer version %s\r\n', version);
            for p = 1:numel(paramNames)
                fd = paramNames{p};
                if strcmp(fd, 'matmot_version'), continue, end
                val = params.(fd);
                if isnumeric(val)
                    fmt = '%u';
                elseif ischar(val)
                    fmt = '%s';
                end
                fprintf(fid, ['%s=' fmt '\r\n'], fd, val);
            end
            fclose(fid);
        end
        
        function editMeta(filename, params)
            import matmot.FormatSpec
            txt = FormatSpec.readMetaText(filename);
            meta = FormatSpec.parseMetaText(txt);
            % Update parsed meta params with new values
            fields = fieldnames(params);
            for f = 1:numel(fields)
                fd = fields{f};
                meta.(fd) = params.(fd);
            end
            matmot.FormatSpec.writeMeta(filename, meta.matmot_version, meta);
        end
        
        function txt = readMetaText(filename)
            % Return parsed metadata for a .mtv file
            [pth, fn, ext] = fileparts(filename);
            if strcmpi(ext, '.mtv')
                ver = matmot.FormatSpec.versionFromFile(filename);
                useMetaFile = ver(2) >= 2;
                if useMetaFile
                    [pth, fn, ext] = fileparts(filename);
                    filename = fullfile(pth, [fn '.meta']);
                end
            elseif strcmpi(ext, '.meta')
                useMetaFile = true;
            end 
            fid = fopen(filename, 'r'); 
            if useMetaFile
                sz = [1, inf];
            else
                sz = [1, matmot.FormatSpec.HEADER_LENGTH];
            end
            txt = fread(fid, sz, '*char');
            fclose(fid);
            txt = strtrim(txt);
        end
        
        function meta = parseMetaText(txt)
            lines = strsplit(txt, '\r\n');
            meta.matmot_version = matmot.FormatSpec.parseVersion(lines{1});
            for i = 2:numel(lines)
                ln = lines{i};
                if ~isempty(ln)
                    idx = strfind(ln, '=');
                    if isempty(idx)
                        warning('matmot:loadMtvFile:invalidMeta', 'Invalid line in metadata file: "%s".', ln);
                    else
                        field = ln(1:idx-1);
                        val = ln(idx+1:end);
                        % Try to parse value as numeric type. If data lost in
                        % conversion then leave as string
                        valNum = str2num(val);
                        if strcmp(num2str(valNum), val)
                            val = valNum;
                        end
                        meta.(field) = val;
                    end
                end
            end
            minorVersion = meta.matmot_version(2);
            if minorVersion < 0.1, meta.n_rigid_bodies = 1; end 
            if minorVersion < 0.2, meta.composite_file = 'false'; end
        end
        
        function v = parseVersion(verTxt)
            strs = strsplit(verTxt);
            verStr = strs{end};
            verElStr = strsplit(verStr, '.');
            assert(numel(verElStr)==3, 'matmot:FormatSpec:noFileVersion', ...
                'Failed to interpret version string "%s"', verStr);
            for i = 1:3
                v(i) = str2num(verElStr{i});
            end
        end
        
        function [v, meta] = versionFromFile(filename)
            import matmot.FormatSpec
            [pth, fn] = fileparts(filename);
            filenameMeta = fullfile(pth, [fn '.meta']);
            fileVerPre02 = ~exist(filenameMeta, 'file');
            if fileVerPre02
                fid = fopen(filename, 'r');
                buff = fread(fid, [1, FormatSpec.HEADER_LENGTH], '*char');
                txt = strtrim(buff);
            else
                fid = fopen(filenameMeta, 'r');
                txt = fread(fid, [1, inf], '*char');
            end
            meta = matmot.FormatSpec.parseMetaText(txt);
            v = meta.matmot_version;
        end
        
        function n = fileNFrames(filename)
            import matmot.FormatSpec
            [v, meta] = FormatSpec.versionFromFile(filename);
            useMetaFile = v(2) > 1;
            if useMetaFile && isfield(meta, 'n_frames')
                n = meta.n_frames;
            else
                if useMetaFile
                    offset = 0;
                else
                    offset = FormatSpec.HEADER_LENGTH;
                end
                fid = fopen(filename, 'r');
                fseek(fid, 0, 'eof');
                dataSize  = ftell(fid) - offset;
                fclose(fid);
                nrb = meta.n_rigid_bodies;
                nm = meta.n_markers;
                frameBytes = FormatSpec.bytesPerFrame(nrb, nm);
                n = dataSize/frameBytes;
                if rem(n, 1) ~= 0
                    error('Data size %u is incorrect for %u rigid bodies and %u markers', dataSize, nrb, nm);
                end
            end
        end
        
    end
    
end

