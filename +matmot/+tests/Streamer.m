classdef Streamer < matlab.unittest.TestCase
    %STREAMER tests, using real or simulated data acquisition
    
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
                self.streamer.deleteFiles();
            end
        end
    end
    
    methods (Test)
        
        function OutputFilesExist(self)
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
        
        function FrameCorrectNBytes(self)
            
        end
        
        function DataWritten(self)
            % Check that the acquired data is faithfully written to file
            
            streamer = self.createStreamer();
            addlistener(streamer, 'postGetFrame', @(src, event) gatherData)
            
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
            
            self.assertTrue(~isempty(fileData), 'Data file was empty');
            self.verifyTrue(isequal(data',fileData), 'File contents did not match acquired bytes');
        end
        
        function DataRead(self)
            % Check that acquired data is faithfully loaded
            streamer = self.createStreamer();
            
            encodings = matmot.Consts.encodings();
            for f = 1:numel(encodings)
                data.(encodings(f).field) = [];
            end
            
            function gatherData()
                % Append current values of data fields as rows
                for f = 1:numel(encodings)
                    field = encodings(f).field;
                    data.(field) = [data.(field); streamer.(field)];
                end
            end
            
            addlistener(streamer, 'postGetFrame', @(src, event) gatherData)
            self.streamer.start();
            pause(1);
            self.streamer.finish();
            filepath = self.streamer.getFiles().data;
            dataLoaded = matmot.loadMtvFile(filepath);
            
            for f = 1:numel(encodings)
                field = encodings(f).field;
                valOriginal = data.(field);
                valLoaded = dataLoaded.(field);
                msg = sprintf('Data field "%s" was empty', field);
                self.assertTrue(~isempty(valOriginal), msg);
                msg = sprintf('Loaded values of field "%s" values different from acquisition values', field);
                self.verifyTrue(isequal(valOriginal, valLoaded), msg)
            end
            

        end
    end
end

