classdef FormatSpec
    %FormatSpec consts and functions defining the matmot file format
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FILE FORMAT DESCRIPTION:
    %
    % The matmot data format is separated into two files:
    % ".meta" - text file containing metadata fields and values
    % ".mtv"  - headerless binary file containing tracking data
    %
    % The MTV file byte structure consists of a series of "records", with
    % the bytes in each record containing data for one camera frame. The
    % size of a record is determined by the number of rigid bodies and
    % markers configured in the Streamer object (this number cannot be
    % changed while recording, so record size does not vary within one .mtv
    % file).
    %
    % The size of a record can be calculated from this formula:
    % NBYTES = 16 + (NUM_RIGID_BODIES * 33) + (NUM_MARKERS * 20)
    %
    %
    % The record comprises a fixed sequence of data fields, each field
    % containing a scalar numerical value with a specific encoding. The
    % fields fall into two categories: "basic" fields which always occur
    % once in each record, and rigid-body or marker fields, which are
    % repeated for each of the rigid bodies and markers that was recorded.
    %
    % FILE RECORD FIELDS:
    %
    % 1) Basic fields: these occur at the start of every record
    % irrespective of how many rigid bodies or markers were configured.
    %
    %   "frameIdx"          (int32)     frame index assigned by Motive
    %   "frameTimestamp"    (double)    frame timestamp assigned by Motive
    %   "frameLatency"      (single)    not currently used
    %
    % 2) Rigid-body fields: this set of fields begins immediately after the
    % basic fields, and repeats once for each rigid body.
    %
    %   "rbx"               (single)    x-position
    %   "rby"               (single)    y-position
    %   "rbz"               (single)    z-position
    %
    %   "rbqx"              (single)    quaternion 1
    %   "rbqy"              (single)    quaternion 2
    %   "rbqz"              (single)    quaternion 3
    %   "rbqw"              (single)    quaternion 4
    %
    %   "rbError"           (single)    error in rigid body solution
    %   "rbTracker"         (uint8)     boolean indicator for whether rigid
    %                                   body was successfully tracked
    %
    %  3) Marker fields: this set of fields begins immediately after the 
    %  rigid-body fields, and repeats once for each marker.
    %   
    %   "mx"                (single)    x-position
    %   "my"                (single)    y-position
    %   "mz"                (single)    z-position
    %
    %   "msz"               (single)    marker size
    %   "mres"              (single)    marker residual
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Static)
        
        % These UPPERCASE methods are really just constants, implemented
        % as functions because matlab doesn't allow importing of class
        % constants.
        
        function val = BASIC_FIELD_NAMES() 
            val = {
            'frameIdx'
            'frameTimestamp'
            'frameLatency'};
        end
        
        function val = BASIC_FIELD_ENCODINGS() 
            val = {
            'int32'
            'double'
            'single'};
        end
        
        function val = RB_FIELD_NAMES()
            val = {
            'rbx'
            'rby'
            'rbz'
            'rbqx'
            'rbqy'
            'rbqz'
            'rbqw'
            'rbError'
            'rbTracked'};
        end
        
        function val = RB_FIELD_ENCODINGS()
            val = {
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'single'
            'uint8'};
        end
        
        function val = MARKER_FIELD_NAMES()
            val = {
            'mx'
            'my'
            'mz'
            'msz'
            'mres'};
        end
        
        function vals = MARKER_FIELD_ENCODINGS() 
            vals = {
            'single'
            'single'
            'single'
            'single'
            'single'};
        end
        
        function val = HEADER_LENGTH()
            val = 2^14;
        end
        
        function val = VERSION()
            val = '0.2.3';
        end
        
        function ver = Version()
            ver = matmot.Version(matmot.FormatSpec.VERSION);
        end
        
        function S = rbFields()
            % Field encodings as a struct array
            import matmot.FormatSpec.*
            [byteInds, nBytes] = getByteInds(RB_FIELD_ENCODINGS);
            S = struct( ...
                'name', RB_FIELD_NAMES, ...
                'byte_inds', num2cell(byteInds'), ...
                'n_bytes', num2cell(nBytes'), ...
                'encoding', RB_FIELD_ENCODINGS);
        end
        
        function S = markerFields()
            % Marker field encodings as a struct array
            import matmot.FormatSpec.*
            [byteInds, nBytes] = getByteInds(MARKER_FIELD_ENCODINGS);
            S = struct( ...
                'name', MARKER_FIELD_NAMES, ...
                'byte_inds', num2cell(byteInds'), ...
                'n_bytes', num2cell(nBytes'), ...
                'encoding', MARKER_FIELD_ENCODINGS);
        end
        
        function S = basicFields()
            % Marker field encodings as a struct array
            import matmot.FormatSpec.*
            [byteInds, nBytes] = getByteInds(BASIC_FIELD_ENCODINGS);
            S = struct( ...
                'name', BASIC_FIELD_NAMES, ...
                'byte_inds', num2cell(byteInds'), ...
                'n_bytes', num2cell(nBytes'), ...
                'encoding', BASIC_FIELD_ENCODINGS);
        end
        
        function nBytes = bytesBasic()
            % Number of bytes for basic frame fields
            fields = matmot.FormatSpec.basicFields();
            nBytes = sum([fields.n_bytes]);
        end
        
        function nBytes = bytesPerRb()
            % Number of bytes needed to encode data of rigid body
            fields = matmot.FormatSpec.rbFields();
            nBytes = sum([fields.n_bytes]);
        end
        
        function nBytes = bytesPerMarker()
            % Number of bytes needed to encode data of one marker
            fields = matmot.FormatSpec.markerFields();
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
            idx0 = 1;
            for n = 1:numel(encodings)
                encoding = encodings{n};
                nBytes(n) = matmot.FormatSpec.encodingNBytes(encoding);
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
            import matmot.FormatSpec.*
            nBytes = ...
                bytesBasic() + ...
                nRbs*bytesPerRb() + ...
                nMarkers*bytesPerMarker();
        end
        
        function writeMeta(filename, version, params)
            paramNames = fieldnames(params);
            fid = fopen(filename, 'w');
            fprintf(fid, 'MatMot Streamer version %s\r\n', version.toString());
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
            import matmot.FormatSpec.*
            
            meta = readMeta(filename);
            % Update parsed meta params with new values
            fields = fieldnames(params);
            for f = 1:numel(fields)
                fd = fields{f};
                meta.(fd) = params.(fd);
            end
            matmot.writeMeta(filename, meta.matmot_version, meta);
        end
        
        function [meta, ver] = parseMetaText(txt)
            import matmot.FormatSpec.*
            import matmot.Version
            lines = strsplit(txt, '\r\n');
            meta.matmot_version = parseVersion(lines{1});
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
            ver = Version(meta.matmot_version);
            if ver < '0.1.0', meta.n_rigid_bodies = 1; end
            if ver < '0.2.0', meta.composite_file = 'false'; end
            if ver < '0.2.2', meta.rbx_inverted = 'true'; end
        end
        
        function meta = readMetaPartFiles(mergedFilename)
            import matmot.FormatSpec.readMeta
            [pth, fn] = fileparts(mergedFilename);
            metaDir = fullfile(pth, [fn '_meta']);
            tmp = dir(fullfile(metaDir, '*.meta'));
            metafiles = {tmp.name};
            nFiles = numel(metafiles);
            txt = {};
            meta = {};
            fileInds = [];
            for f = 1:nFiles
                fn = fullfile(metaDir, metafiles{f});
                meta{f, 1} = readMeta(fn);
                fileInds(f) = sscanf(metafiles{f}, '%u.meta');
            end
            [~, iSort] = sort(fileInds);
            txt = txt(iSort);
            meta = meta(iSort);
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
        
        function [meta, ver, offset] = readMeta(filename)
            import matmot.FormatSpec.*
            [pth, fn] = fileparts(filename);
            filenameMeta = fullfile(pth, [fn '.meta']);
            useMetaFile = exist(filenameMeta, 'file');
            if useMetaFile
                fid = fopen(filenameMeta, 'r');
                metaTxt = fread(fid, [1, inf], '*char');
                offset = 0;
            else
                fid = fopen(filename, 'r');
                buff = fread(fid, [1, HEADER_LENGTH], '*char');
                metaTxt = strtrim(buff);
                offset = HEADER_LENGTH;
            end
            fclose(fid);
            [meta, ver] = parseMetaText(metaTxt);
            
            fid = fopen(filename, 'r');
            fseek(fid, 0, 'eof');
            dataSize  = ftell(fid) - offset;
            fclose(fid);
            nrb = meta.n_rigid_bodies;
            nm = meta.n_markers;
            frameBytes = bytesPerFrame(nrb, nm);
            nFrames = dataSize/frameBytes;
            if rem(nFrames, 1) ~= 0
                warning('matmot:badFileSize', 'Data size %u is incorrect for %u rigid bodies and %u markers', dataSize, nrb, nm);
            end
            nFrames = floor(nFrames);
            
            if isfield(meta, 'n_frames') && meta.n_frames ~= nFrames
                warning('matmot:badMetaNFrames', ...
                    'Frame count mismatch: file contains %u frames but meta.n_frames=%u.', ...
                    nFrames, meta.n_frames)
            end
            meta.n_frames = nFrames;
        end
        
    end
    
end

