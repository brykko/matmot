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
        
    end
    
end

