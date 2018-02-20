classdef MotiveStreamer < matlab.unittest.TestCase
    %MOTIVESTREAMER tests, using simulated data acquisition
    
    properties
        streamer
    end
    
    methods
        function createStreamer(self, varargin)
            pth = fileparts(mfilename('fullpath'));
            tempFilePath = fullfile(pth, 'motive_test');
            self.streamer = matmot.MotiveStreamer( ...
                'simulate', true, ...
                'fileName', tempFilePath, ...
                varargin{:});
        end
    end
    
    methods (TestMethodTeardown)
        function deleteTempFiles(self)
            self.streamer.deleteFiles();
        end
    end
    
    methods (Test)
        
        function OutputFilesExist(self)
            % Check the three output files exist with the names they're
            % supposed to have
            self.createStreamer();
            self.streamer.start();
            self.streamer.finish();
            
            filepaths = self.streamer.getFiles();
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
        
        function DataWritten(self)
            % Check that the acquired data is faithfully written to file
            
            self.createStreamer()
            streamer = self.streamer;
            addlistener(streamer, 'postGetFrame', @(src, event) gatherData)
            
            % Callback appends frame data to array every time a new frame
            % is captured
            data = zeros(0, 'uint8');
            function gatherData()
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
        
    end
    
end

