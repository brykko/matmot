classdef Consts < handle
    %CONSTS file format-related constants and functions
    
    properties (Constant)
        
        FIELDS = {
            'frameIdx'
            'frameTimestamp'
            'frameLatency'
            'pos'
            'rot'
            'posError'
            'posTracked'};
        
        ENCODINGS = {
            'int32'
            'double'
            'single'
            '3*single'
            '4*single'
            'single'
            'uint8'};
        
        MARKER_FIELDS = {
            'mx'
            'my'
            'mz'
            'msz'
            'mres'};
        
        MARKER_ENCODINGS = {
            'single'
            'single'
            'single'
            'single'
            'single'
            }
        
    end
    
    methods (Static)
        
        function S = encodings()
            % Field encodings as a struct array
            import matmot.Consts
            [byteInds, nBytes] = Consts.getByteInds(Consts.ENCODINGS);
            S = struct( ...
                'field', Consts.FIELDS, ...
                'byte_inds', num2cell(byteInds'), ...
                'n_bytes', num2cell(nBytes'), ...
                'encoding', Consts.ENCODINGS);
        end
        
        function S = markerEncodings()
            % Marker field encodings as a struct array
            import matmot.Consts
            [byteInds, nBytes] = Consts.getByteInds(Consts.MARKER_ENCODINGS);
            S = struct( ...
                'field', Consts.MARKER_FIELDS, ...
                'byte_inds', num2cell(byteInds'), ...
                'n_bytes', num2cell(nBytes'), ...
                'encoding', Consts.MARKER_ENCODINGS);
        end
        
        function nBytes = bytesPerMarker()
            % Number of bytes needed to encode all of a marker's fields
            import matmot.Consts
            encodings = Consts.markerEncodings();
            nBytes = sum([encodings.n_bytes]);
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
            [encoding, n] = matmot.Consts.convertEncoding(encoding);
            x = zeros(1, n, encoding);
            nBytes = numel(typecast(x, 'uint8'));
        end
        
        function [inds, nBytes] = getByteInds(encodings)
            % Byte indices and numbers 
            import matmot.Consts
            idx0 = 1;
            for n = 1:numel(encodings)
                encoding = encodings{n};
                nBytes(n) = Consts.encodingNBytes(encoding);
                inds(n) = idx0;
                idx0 = idx0 + nBytes(n);
            end
        end
        
        function nBytes = bytesPerFrame(nMarkers)
            %BYTESPERFRAME number of bytes needed to encode one frame
            %
            % NBYTES = BYTESPERFRAME(NMARKERS) returns the number of bytes 
            % NBYTES needed to encode a frame containing the maximum 
            % marker count NMARKERS.
            import matmot.Consts
            encodings = Consts.encodings();
            nBytesBasic = sum([encodings.n_bytes]);
            nBytes = nBytesBasic + nMarkers*Consts.bytesPerMarker();
        end
        
    end
    
end

