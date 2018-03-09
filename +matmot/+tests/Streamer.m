classdef Streamer < matlab.unittest.TestCase
    %STREAMER tests, using real or simulated data acquisition
    
    %TODO add tests for these:
    % - check CustomDisplay props all exist
    
    properties (Constant)
        % Enable this parameter to use simulated (randomly generated) data.
        % If disabled, the Streamer will attempt to connect to the NatNet
        % server at localhost.
        USE_SIMULATION = false;
    end
    
    properties
        streamer
    end
    
    methods
        function streamer = createStreamer(self, varargin)
            pth = fileparts(mfilename('fullpath'));
            tempFilePath = fullfile(pth, 'streamer_test');
            streamer = matmot.Streamer( ...
                'simulate', self.USE_SIMULATION, ...
                'fileName', tempFilePath, ...
                varargin{:});
            self.streamer = streamer;
        end
    end
    
    methods (TestMethodTeardown)
        function deleteTempFiles(self)
            if ~isempty(self.streamer) && isvalid(self.streamer)
                self.streamer.finish();
                self.streamer.deleteFiles();
            end
        end
    end
    
    methods (Test)
        
        function outputFilesExist(self)
            % Check the three output files exist with the names they're
            % supposed to have
            streamer = self.createStreamer();
            streamer.start();
            streamer.finish();
            
            filepaths = streamer.getFiles();
            fields = fieldnames(filepaths);
            self.assertTrue( ...
                isstruct(filepaths) && numel(fields)==3, ...
                'Expected filepaths structure with 3 fields');
            
            for f = 1:numel(fields)
                fd = fields{f};
                pth = filepaths.(fd);
                fileExists = exist(pth, 'file') > 0;
                msg = sprintf('Failed to find %s file "%s"', fd, pth);
                self.verifyTrue(fileExists, msg);
            end
        end
        
        function frameCorrectNBytes(self)
            % Check that the number of encoded bytes from a frame match the
            % expected number
            nRbs = [0, 1, 2];
            nMarkers = [10, 0, 10];
            for n = 1:3
                streamer = self.createStreamer( ...
                    'nRigidBodies', nRbs(n), ...
                    'nMarkers', nMarkers(n));
                streamer.start();
                pause(0.1);
                self.assertFalse(isempty(streamer.firstFrameIdx), ...
                    'Streamer has no frame data yet')
                bytes = streamer.currentFrameToBytes();
                streamer.finish();
                self.verifyTrue(numel(bytes)==streamer.nBytesPerFrame, ...
                    'Mismatch between expected and actual number of frame bytes');
            end
        end
        
        function frameCorrectEncodings(self)
            % Check that all fields in a frame are of the expected encoding
            import matmot.FormatSpec
            streamer = self.createStreamer();
            
            streamer.start();
            checkAllFields();
            pause(0.1);
            checkAllFields();
            
            self.assertFalse(isempty(streamer.firstFrameIdx), ...
                'Streamer has no frame data yet')
            
            function checkAllFields()
                import matmot.FormatSpec;
                fields = [
                    FormatSpec.basicFields()
                    FormatSpec.rbFields()
                    FormatSpec.markerFields()];
                
                nRows = [
                    ones(size(FormatSpec.basicFields()))
                    streamer.nRigidBodies * ones(size(FormatSpec.rbFields()))
                    streamer.nMarkers * ones(size(FormatSpec.markerFields()))
                    ];
                
                for f = 1:numel(fields)
                    name = fields(f).name;
                    [encoding, nCols] = FormatSpec.convertEncoding(fields(f).encoding);                 
                    val = streamer.(name);
                    
                    % Check that class matches the required encoding
                    matchClass = isa(val, encoding);
                    msg = sprintf('Expected field "%s" to contain a value of class %s, but was %s instead', ...
                        name, encoding, class(val));
                    self.verifyTrue(matchClass, msg);
                    
                    % Check the number of columns matches the encoding
                    matchSize = isequal(size(val), [nRows(f), nCols]);
                    msg = sprintf('Expected field "%s" to contain array of size [%u, %u], but was [%u, %u] instead', ...
                        name, nRows(f), nCols, size(val, 1), size(val, 2));
                    self.verifyTrue(matchSize, msg);
                end
            end
            
        end
        
        function dataWritten(self)
            % Check that the acquired data is faithfully written to file
            
            streamer = self.createStreamer();
            addlistener(streamer, 'postGetFrame', @(src, event) gatherData);
            
            % Callback appends frame data to array every time a new frame
            % is captured
            data = zeros(0, 'uint8');
            function gatherData()
                % Append bytes as row vector
                data = [data streamer.currentFrameToBytes()];
            end
            
            % Capture some data
            streamer.start();
            pause(1);
            streamer.finish();
            
            % Check the data file
            filepath = streamer.getFiles().data;
            fid = fopen(filepath, 'r');
            fread(fid, streamer.HEADER_LENGTH);
            fileData = fread(fid, '*uint8');
            fclose(fid);
            
            self.assertFalse(isempty(fileData), 'Data file was empty');
            self.verifyTrue(isequal(data',fileData), 'File contents did not match acquired bytes');
        end
        
        function dataRead(self)
            
            fields = [
                matmot.FormatSpec.basicFields()
                matmot.FormatSpec.rbFields()];
            
            function gatherData()
                % Append current values of data fields as rows
                for f = 1:numel(fields)
                    name = fields(f).name;
                    data.(name) = [data.(name); self.streamer.(name)];
                end
            end
            
            nMarkers = [0, 20]; % test both with and without markers
            
            for m = 1:2
                % Check that acquired data is faithfully loaded
                streamer = self.createStreamer('nMarkers', nMarkers(m));
                for f = 1:numel(fields)
                    data.(fields(f).name) = [];
                end
                addlistener(streamer, 'postGetFrame', @(src, event) gatherData);
                
                self.streamer.start();
                pause(1);
                self.streamer.finish();
                
                filepath = self.streamer.getFiles().data;
                dataLoaded = matmot.loadMtvFile(filepath);
                for f = 1:numel(fields)
                    name = fields(f).name;
                    valOriginal = data.(name);
                    valLoaded = dataLoaded.(name);
                    msg = sprintf('Data field "%s" was empty', name);
                    self.assertFalse(isempty(valOriginal), msg);
                    msg = sprintf('Loaded values of field "%s" values different from acquisition values', name);
                    self.verifyTrue(isequal(valOriginal, valLoaded), msg)
                end
            end
        end
        
        function pauseStopsAqcuisition(self)
            streamer = self.createStreamer();
            streamer.start();
            pause(0.1);
            self.assertTrue(streamer.nFramesAcquired > 0, 'Streamer failed to capture any frames')
            
            % Pause the streamer, wait a while and check that no more
            % frames have been acquired
            streamer.pause();
            n1 = streamer.nFramesAcquired;
            pause(0.1);
            n2 = streamer.nFramesAcquired;
            msg = sprintf('Number of acquired frames increased from %u to %u while paused', n1, n2);
            self.verifyTrue(n1==n2, msg);
        end
        
        function resumeContinuesAcquisision(self)
            streamer = self.createStreamer();
            streamer.start();
            streamer.pause();
            n1 = streamer.nFramesAcquired;
            streamer.start();
            pause(0.1);
            n2 = streamer.nFramesAcquired;
            self.verifyTrue(n1~=n2, 'Number of frames did not change after resuming')
        end
        
        function loadSampleData(self)
            rootPath = fileparts(fileparts(mfilename('fullpath')));
            sampleFilePath = fullfile(rootPath, '+samples', 'motive_stream.mtv');
            matmot.loadMtvFile(sampleFilePath);
        end
        
        function noDataTimeout(self)
            function callback(~, ~)
                n = n+1;
            end
            
            % With long timeout, check the event isn't notified
            n = 0;
            streamer = self.createStreamer();
            streamer.noDataTimeout = 1;
            addlistener(streamer, 'noData', @callback);
            streamer.start();
            pause(0.1);
            streamer.finish();
            self.verifyTrue(n==0, 'Timeout warning was notified');
            
            % With short timeout, it should be notified multiple times
            n = 0;
            streamer = self.createStreamer();
            streamer.noDataTimeout = 1/300;
            addlistener(streamer, 'noData', @callback);
            streamer.start();
            pause(0.1);
            streamer.finish();
            self.verifyTrue(n>1, 'Timeout warning was not notified');
            
            
        end
        
    end
end