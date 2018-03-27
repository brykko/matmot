classdef Streamer < handle & matlab.mixin.CustomDisplay
    %STREAMER stream OptiTrack Motive rigid body and marker data to disk.
    %
    % Streamer objects stream data from Motive, saving the received data 
    % frames to a binary file.
    %
    % S = STREAMER() creates a new Streamer object S with default parameter
    % settings.
    %
    % S = STREAMER(PRM, VAL, ... ) creates a new Streamer object S and
    % initializes it with the specified parameter-value pair arguments. The 
    % available parameters are as follows:
    %
    %   'hostIp' (default '127.0.0.1') host IP address
    %
    %   'fileName' (default 'motive_stream.mtv') name of output file
    %
    %   'frameIncrement' (default 1) save received frames with this
    %       interval. A value of 1 saves every frame; a value of 3 saves
    %       every third frame, etc.
    %
    %   'nRigidBodies' (default 1) number of rigid bodies to record in
    %   file. The x/y/z position coordinates and quaternions are saved for
    %   each rigid body. Rigid bodies are saved in the same order as they
    %   are configured in Motive.
    %
    %   'nMarkers' (default 0) maximum number of labelled markers to record 
    %   in file. Markers are recorded independently of rigid bodies, so
    %   this option can be useful if you have markers of interest that are 
    %   not members of a rigid body (e.g. IR-LEDs for syncing). The 
    %   x/y/z coordinates of markers received in frames from the NatNet
    %   client will be written to the .mtv file, up to the limit specified
    %   by nMarkers. Increasing the value of nMarkers will reduce the
    %   probability that markers of interest are not recorded, but will 
    %   increase the processing load and may result in more dropped frames.
    %
    %   'writeToFile' (default TRUE) specifies whether acquired frames will
    %   be written to the .mtv file or not. If the value is FALSE, the
    %   Streamer will still acquired data and notify events as normal, but
    %   no data will be written to file.
    %
    %   'noDataTimeout' (default 1.0) the number of seconds to wait to
    %   receive new frames from Motive before displaying a warning and
    %   notifying the 'noData' event.
    %
    %   'untrackedWarnTreshold' (default 30) the minimum percentage of
    %   frames for which a rigid body may be untracked before issuing a
    %   warning.
    %
    %   'untrackedWarnInterval' (default 10) the time interval in seconds 
    %   for checking how many untracked frames each rigid body has 
    %   accumulated.
    %
    % --------------------------------------------------------------------
    % STREAMING
    %
    % S.START() begins data acquisition.
    %
    % S.PAUSE() temporarily stops acquisition. Call start() again to
    % resume.
    %
    % S.FINISH() finishes acquisition and closes the data file.
    
    % TODO: replace overly verbose loss-of-tracking warnings with something
    % more useful
    
    properties (SetAccess = protected, Transient)
        NNClient
    end
    
    properties
        hostIP              (1,:) char = '127.0.0.1'
        
        % Ouput file props
        fileName            (1,:) char = 'motive_stream.mtv'
        
        % Aqcuisition settings
        writeToFile         (1,1) logical = true
        frameIncrement      (1,1) double {mustBeInteger, mustBePositive} = 1
        writeBufferNFrames  (1,1) double {mustBeInteger, mustBePositive} = 120
        nRigidBodies        (1,1) double {mustBeInteger, mustBeNonnegative} = 1
        nMarkers            (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        noDataTimeout       (1,1) double {mustBeNonnegative} = 1;
        untrackedWarnTreshold (1,1) double {mustBeNonnegative} = 30;
        untrackedWarnInterval (1,1) double {mustBeNonnegative} = 10;
        
        % Misc
        debug               (1,1) logical = false
        simulate            (1,1) logical = false
    end
    
    properties (SetAccess = protected)
        % Aqcuisition results
        frameRate
        nFramesAcquired = 0;
        nFramesDropped = 0;
        nFramesInBuffer = 0;
        
        meanGetFrameTime
        firstFrame = true;
        firstFrameIdx
        firstFrameTimestamp = 0;
        frameIdx = int32(0);
        lastTrackedFrameIdx
        
        nUntrackedFramesSinceWarn
        nTotalFramesSinceWarn = 0;
        
        % Position variables: each contains the value from the most recent
        % frame
        rbx
        rby
        rbz
        rbError
        rbTracked
        
        % Quaternion angles
        rbqx
        rbqy
        rbqz
        rbqw
        
        % Marker info (don't initialize here; dims are unknown)
        mx
        my
        mz
        msz
        mres

        frameTimestamp = 0;
        frameLatency = single(0)
    end
    
    properties (SetAccess = protected, Dependent)
        timeElapsed
        pos
        rot
    end
    
    properties (Constant)
        HEADER_LENGTH = 2^14
        VERSION = '0.1.0'
    end
    
    properties (SetAccess = protected)
        % Streaming state
        clientInitialized = false;
        streamingInitialized = false;
        streaming = false;
        streamingFinished = false;
        
        % Data file
        fileOpen = false;
        fid
        dataDir
        writeBuffer = uint8.empty()
        writeBufferTmp
        logger
        nWriteCycles = 0;
        meanWriteTime
    end
    
    properties (SetAccess = protected, Transient)
        frameReadyListener
        noDataTimer
        rbUntrackedTimer
    end
    
    events
        % EVENTS
        %
        % Attach listeners to these salient events to implement custom
        % callbacks
        
        preGetFrame      % Immediately before a new frame is written
        postGetFrame     % Immediately after a new frame is written
        droppedFrame     % One or more dropped frames are detected
        duplicateFrame   % The registered frame is identical to the last
        untrackedFrame   % An untracked frame is registered
        trackingLost     % A period of untracked frames begins
        trackingResumed  % A period of untracked frames ends
        noData           % No frames received during specified timeout period
    end
    
    methods
        
        function self = Streamer(varargin)
            % STREAMER constructor
            %
            % S = STREAMER(DLLPATH) creates a new Streamer
            % object S using the NatNet DLL at the path specified by
            % DLLPATH.
            %
            % S = STREAMER(DLLPATH, PRM, VAL, ... ) creates a new
            % Streamer object S and initializes it with specified
            % parameter/values pairs. See main docstring for available
            % parameters.
            self.parseInputs(varargin{:});
            self.initializeLogging();
            self.initializeClient();
            self.logger.i('Initialization complete. Call start() to begin streaming.');
        end
        
        function start(self)
            % START - begin or resume streaming
            %
            % S.START() starts acquistion of data to disk. Acquisition
            % continues until the pause() or finish() methods are called.
            %
            % Calling S.START() after previously calling pause() will
            % resume streaming.
            
            % Check finish() hasn't been called
            if self.streamingFinished
                self.logger.w('Cannot resume streaming after calling finish()');
                return;
            end
            
            % Initialize all of the vars to the appropriate starting values
            if ~self.streamingInitialized
                self.initializeStreaming();
                self.logger.i('Starting streaming. Call finish() to stop.');
                if self.simulate
                    % Simulate frame ready events with a 120 Hz timer
                    tim = timer( ...
                        'Period', 1/120, ...
                        'TimerFcn', @(~,~) callback(self), ...
                        'ExecutionMode', 'fixedRate');
                    start(tim);
                    self.frameReadyListener = tim;
                else
                    % Attach a listener to execute whenever the NatNet client
                    % notifies the event "OnFrameReady2"
                    self.frameReadyListener = addlistener(self.NNClient, 'OnFrameReady2', @(~,~) callback(self));
                end
            end
            
            start(self.noDataTimer)
            start(self.rbUntrackedTimer)
            self.streaming = true;
            
            function callback(streamer)
                if streamer.streaming
                    streamer.getFrame();
                end
            end
            
        end
        
        function pause(self)
            % PAUSE - temporarily pause streaming
            %
            % S.PAUSE() stops Streamer object S from receiving new tracking
            % data.
            if ~self.streaming
                self.logger.w('Not currently streaming. Calling pause() will have no effect!');
            end
            stop(self.noDataTimer);
            stop(self.rbUntrackedTimer);
            self.streaming = false;
        end
        
        function finish(self)
            % FINISH - finish streaming
            %
            % S.FINISH() stops acquiring data, closes the data file and
            % saves Streamer object S to a .mat file.
            
            % Some checks:
            % finish() may be called before streaming has begun
            if ~self.streamingInitialized
                self.logger.w('Streaming has not been initialized. Calling finish() will have no effect!');
                return;
            end
            
            stop(self.noDataTimer);
            stop(self.rbUntrackedTimer);
            
            % finish() may be called after a previous call to finish()
            if self.streamingFinished
                self.logger.w('Streaming has already finished. Calling finish() will have no effect!');
                return;
            end
            
            % For callback mode, delete the FrameReady listener
            if ~isempty(self.frameReadyListener) && isvalid(self.frameReadyListener)
                delete(self.frameReadyListener);
            end
            
            % Flush write buffer
            if ~isempty(self.writeBuffer)
                self.logger.d('%u frames still in buffer at end of streaming. Writing now...', self.nFramesInBuffer);
                self.flushBuffer();
            end
            
            self.closeOutputFile();
            
            % Done! Generate summary and clean up
            self.logger.i('Streaming finished. Total number of dropped frames: %u. Mean getFrame time = %.2f ms, flushBuffer time = %.2f ms.', ...
                self.nFramesDropped, self.meanGetFrameTime*1e3, self.meanWriteTime*1e3);
            
            % Save self to .mat file
            filePath = self.getFiles().streamer;
            streamerObj = self;
            save(filePath, 'streamerObj');
            self.logger.i('Saved Streamer obj in file %s', filePath);
            
            self.streaming = false;
            self.streamingFinished = true;
            self.closeLogFile();
        end
        
        function paths = getFiles(self)
            % GETFILES return filepaths associated with Streamer instance
            [~, baseName, ~] = fileparts(self.fileName);
            paths = struct( ...
                'data', fullfile(self.dataDir, self.fileName), ...
                'log', fullfile(self.dataDir, [baseName '.log']), ...
                'streamer', fullfile(self.dataDir, [baseName '_streamer.mat']));
        end
        
        function deleteFiles(self)
            % DELETEFILES remove all files associated with Streamer instance
            
            if self.streaming
                error('Please finish streaming before deleting files')
            end
            
            self.closeOutputFile();
            self.closeLogFile();
            filePaths = self.getFiles();
            fields = fieldnames(filePaths);
            for f = 1:numel(fields)
                fd = fields{f};
                pth = filePaths.(fd);
                fprintf('Deleting %s file "%s"\n', fd, pth);
                if exist(pth, 'file')
                    delete(pth)
                else
                    warning('Could not delete %s file "%s": file not found.', fd, pth)
                end
            end
        end
        
        function val = get.timeElapsed(self)
            val = self.frameTimestamp - self.firstFrameTimestamp;
        end
        
    end
    
    methods (Hidden)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INTERNAL METHODS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function initializeLogging(self)
            
            import logging.Level
            
            % Get and reset the class logger
            logger = logging.getLogger(class(self));
            logger.reset();
            logger.useParentHandlers = false;
            logger.level = Level.DEBUG;
            
            % Set up the handlers
            handler = logging.ConsoleHandler();
            handler.useColors = false;
            handler.level = Level.INFO;
            logger.addHandler(handler);
            
            filePath = self.getFiles().log;
            handler = logging.FileHandler(filePath, false);
            handler.level = Level.DEBUG;
            handler.logger.level = Level.DEBUG;
            logger.addHandler(handler);
            
            lg = @(varargin) logger.i(varargin{:});
            lg('New Streamer session, date %s, user "%s", PC "%s"', ...
                datestr(now, 'dd/mm/yyyy'), getenv('username'), getenv('computername'));
            
            self.logger = logger;
            logger.v('Logging intitialization complete');
        end
        
        function initializeClient(self)
            % INITIALIZECLIENT - set up the NatNet Client
            %
            % Get the NatNet client.
            % If there's a preexisting client, it means the last call to this function
            % did not complete, so we need to uninitialize it to prevent MATLAB crashing.
            
            self.logger.d('Initializing NatNet client...');
            
            if isempty(self.NNClient)
                % No client currently exists: create a new one
                % Find DLL path
                basepath = fileparts(mfilename('fullpath'));
                dllPath = fullfile(basepath, '+external', 'NatNetML.dll');
                try
                    [~] = NET.addAssembly(dllPath);
                catch err
                    msg = sprintf('Bad DLL path %s (message : "%s")', dllPath, err.message);
                    self.logger.f('%s', msg);
                    error(msg);
                end
                client = NatNetML.NatNetClientML(0); % Input = iConnectionType: 0 = Multicast, 1 = Unicast
                version = client.NatNetVersion();
                self.logger.d('Created new NatNet Client, Version : %d.%d.%d.%d', version(1), version(2), version(3), version(4) );
                self.NNClient = client;
            else
                % Client already exists: uninitialize it before
                % reinitializing
                client = self.NNClient;
                self.logger.d('NatNet Client already exists - uninitializing');
                client.Uninitialize();
                self.clientInitialized = false;
            end
            
            % Connect to a local stream.
            ip = self.hostIP;
            self.logger.d('Initializing NatNet client at IP %s', ip);
            
            % If in simulation mode, don't try to connect to server
            if self.simulate
                self.frameRate = 120;
            else
                result = client.Initialize(ip, ip); % Flg = returnCode: 0 = Success
                
                if result == 0
                    self.logger.d('NatNet initialization successful!')
                else
                    msg = 'NatNet initialization Failed';
                    self.logger.f(msg);
                    error(msg);
                end
                
                % Display a summary of the tracked rigid body
                self.printDescriptions(client.GetDataDescriptions());
                
                % Test - send command/request to Motive
                [bytes, retCode] = client.SendMessageAndWait('FrameRate');
                if(retCode ==0)
                    self.frameRate = typecast(uint8(bytes), 'single');
                end
            end
            
            self.clientInitialized = true;
            
        end
        
        function initializeStreaming(self)
            % INITIALIZESTREAMING - called just before beginning streaming
            
            self.logger.i('Initializing streaming');
            
            % Open output file
            filePath = self.getFiles().data;
            self.logger.v('Opening output data file %s', filePath);
            [self.fid, message] = fopen(filePath, 'w');
            if self.fid == -1
                msg = sprintf('Failed to open output file %s, message "%s"', filePath, message);
                self.logger.f('%s', msg);
                error(msg);
            end
            self.fileOpen = true;
            self.writeHeader();
            
            % Initialize various streaming variables
            self.writeBufferTmp = zeros(1, self.nBytesPerFrame(), 'uint8');
            
            markerFields = {'mx', 'my', 'mz', 'msz', 'mres'};
            for f = 1:numel(markerFields)
                self.(markerFields{f}) = zeros(self.nMarkers, 1, 'single');
            end
            
            rbFields = {'rbx', 'rby', 'rbz', 'rbqx', 'rbqy', 'rbqz', 'rbqw', ...
                'rbError'};
            for f = 1:numel(rbFields)
                self.(rbFields{f}) = zeros(self.nRigidBodies, 1, 'single');
            end
            self.rbTracked = zeros(self.nRigidBodies, 1, 'uint8');
            
            self.lastTrackedFrameIdx = nan(self.nRigidBodies, 1);
            self.nUntrackedFramesSinceWarn = zeros(self.nRigidBodies, 1);
            
            % Start a timer that will emit warnings at regular intervals if
            % no frames are received within a user-defined timeout period.
            % Every call to getFrame() resets this timer, so it will never
            % fire its callback unless getFrame() fails to be called for
            % the whole timeout period
            self.noDataTimer = timer( ...
                'ExecutionMode', 'fixedSpacing', ...
                'Period', self.noDataTimeout, ...
                'StartDelay', self.noDataTimeout, ...
                'TimerFcn', @(~,~) noDataWarning(self));
            
            self.rbUntrackedTimer = timer( ...
                'ExecutionMode', 'fixedSpacing', ...
                'Period', self.untrackedWarnInterval, ...
                'StartDelay', self.untrackedWarnInterval, ...
                'TimerFcn', @(~,~) rbLostWarning(self));
            
            self.streamingInitialized = true;
            
            function noDataWarning(streamer)
                streamer.logger.w('WARNING! No frame data received for %.1g s', streamer.noDataTimeout);
                streamer.notify('noData');
            end
            
            function rbLostWarning(streamer)
                nUntracked = streamer.nUntrackedFramesSinceWarn;
                nTotal = streamer.nTotalFramesSinceWarn;
                for r = 1:self.nRigidBodies
                    prcUntracked = nUntracked(r) / nTotal * 100;
                    if prcUntracked >= self.untrackedWarnTreshold
                        self.logger.w('Rigid body #%u untracked in %.0f%% of frames during last %.1f seconds', ...
                            r, prcUntracked, self.untrackedWarnInterval);
                    end
                    streamer.nUntrackedFramesSinceWarn(r) = 0;
                end
                streamer.nTotalFramesSinceWarn = 0;
            end
            
        end
        
        function delete(self)
            % Destructor: handle any remaining cleanup processes
            self.logger.d('Calling class destructor');
            if self.streamingInitialized && ~self.streamingFinished
                self.finish();
            end
            delete(self.noDataTimer);
            self.closeClient();
            self.closeOutputFile();
            self.closeLogFile();
        end
        
        function writeHeader(self)
            
            function printprm(prm, val, fmt)
                if nargin < 3 || isempty(fmt), fmt = '%s'; end
                assert(~any(isspace(prm)))
                fprintf(self.fid, ['%s=' fmt '\r\n'], prm, val);
            end
            
            % Write parameter settings to header
            fprintf(self.fid, 'MatMot Streamer version %s\r\n', self.VERSION);
            tmp = self.NNClient.NatNetVersion();
            printprm('natnet_version', tmp(1), '%u')
            printprm('matlab_version', version())
            printprm('time_started', datestr(now()))
            printprm('host_ip', self.hostIP)
            printprm('data_dir', self.dataDir)
            printprm('file_name', self.fileName)
            printprm('frame_increment', self.frameIncrement, '%u')
            printprm('writebuff_nframes', self.writeBufferNFrames, '%u')
            printprm('n_rigid_bodies', self.nRigidBodies, '%u')
            printprm('n_markers', self.nMarkers, '%u')
            printprm('simulate', self.simulate, '%u')
            
            % Fill remainder of 16 kb header allocation with white space
            nBytesWritten = ftell(self.fid);
            bytes = repmat(char(32), 1, self.HEADER_LENGTH-nBytesWritten);
            fwrite(self.fid, bytes);
            
        end
        
        function closeOutputFile(self)
            if self.fileOpen
                filePath = fullfile(self.dataDir, self.fileName);
                self.logger.d('Closing file "%s"', filePath);
                fclose(self.fid);
                self.fileOpen = false;
            end
        end
        
        function closeLogFile(self)
            fileHandler = self.logger.getHandlers('logging.FileHandler');
            if ~isempty(fileHandler)
                fileHandler = fileHandler{1};
                fileHandler.close();
                self.logger.removeHandler(fileHandler);
            end
        end
        
        function closeClient(self)
            if ~isempty(self.NNClient)
                if self.clientInitialized
                    self.logger.d('Uninitializing NatNet client');
                    self.NNClient.Uninitialize();
                    self.clientInitialized = false;
                end
                self.NNClient.Dispose();
            end
        end
        
        function n = nBytesPerFrame(self)
            n = matmot.FormatSpec.bytesPerFrame(self.nRigidBodies, self.nMarkers);
        end
        
        function parseInputs(self, varargin)
            
            inp = inputParser();
            inp.CaseSensitive = false;
            inp.PartialMatching = false;
            inp.KeepUnmatched = false;
            
            prms = {
                'fileName'
                'writeToFile'
                'frameIncrement'
                'hostIP'
                'writeBufferNFrames'
                'nMarkers'
                'nRigidBodies'
                'simulate'};
            
            for p = 1:numel(prms)
                prm = prms{p};
                inp.addParameter(prm, self.(prm));
            end
            inp.parse(varargin{:});
            P = inp.Results;
            
            for p = 1:numel(prms)
                self.(prms{p}) = P.(prms{p});
            end
            
            % Determine data directory
            [pth, fn, ext]  = fileparts(self.fileName);
            if isempty(pth), pth = pwd(); end
            self.dataDir = pth;
            self.fileName = [fn '.mtv'];
            
        end
        
        function flushBuffer(self)
            % FLUSHBUFFER write contents of buffer to disk
            if ~isempty(self.writeBuffer)
                self.logger.v('Flushing write buffer');
                try
                    tic()
                    fwrite(self.fid, self.writeBuffer, '*uint8');
                    self.writeBuffer = [];
                    self.nWriteCycles = self.nWriteCycles+1;
                    self.meanWriteTime = runningMean( ...
                        toc(), self.nWriteCycles, self.meanWriteTime);
                catch e
                    warning('Error while trying to write %u frames: "%s"', ...
                        self.nFramesInBuffer, e.message);
                end
            end
        end
        
        function getFrame(self)
            % GETFRAME acquire a single frame.
            %
            % This may be called either within a polling loop, or by a
            % callback triggered on "FrameReady" events
            
            import matmot.RbEventData
            
            tic();
            self.logger.v('getFrame()');
            
            % Restart the no-frame-data warning timer
            stop(self.noDataTimer);
            start(self.noDataTimer);
            
            % Get new frame of data
            if self.simulate
                newFrameIdx = int32(self.nFramesAcquired+1);
                newFrameTime = (self.nFramesAcquired+1) / self.frameRate;
            else
                frame = self.NNClient.GetLastFrameOfData();
                newFrameIdx = frame.iFrame;
                newFrameTime = frame.fTimestamp;
            end
            
            newFrame = true;
            
            if ~self.firstFrame
                % Check for dropped frames
                frameInc = newFrameIdx - self.frameIdx;
                if(frameInc > self.frameIncrement) && ~self.firstFrame
                    nDropped = frameInc - self.frameIncrement;
                    self.nFramesDropped = self.nFramesDropped + nDropped;
                    self.logger.w('Dropped Frame(s) : %d\tLastIdx : %u\tThisIdx : %u', ...
                        nDropped, self.frameIdx, newFrameIdx);
                    self.logger.v('Notifying "droppedFrame"');
                    self.notify('droppedFrame');
                    
                    % Check for duplicate frame
                elseif(frameInc < self.frameIncrement)
                    newFrame = false;
                    self.logger.v('Notifying "duplicateFrame"');
                    self.notify('duplicateFrame');
                    self.logger.v('Duplicate frame detected, ID = ', newFrameIdx);
                end
            end
            
            if newFrame
                
                self.logger.v('Notifying "pregetFrame"');
                self.notify('preGetFrame');
                self.logger.v('New frame detected, ID = ', newFrameIdx);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % UPDATE TRACKING STATE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Record time and index of the first frame
                if self.firstFrame
                    self.firstFrameTimestamp = newFrameTime;
                    self.firstFrameIdx = newFrameIdx;
                end
                
                if self.simulate
                    % Position and angles are random walks
                    self.rbx = self.rbx + randn(self.nRigidBodies, 1, 'single');
                    self.rby = self.rby + randn(self.nRigidBodies, 1, 'single');
                    self.rbz = self.rbz + randn(self.nRigidBodies, 1, 'single');
                    self.rbqx = mod(self.rbqx + randn(self.nRigidBodies, 1, 'single'), 2*pi);
                    self.rbqy = mod(self.rbqy + randn(self.nRigidBodies, 1, 'single'), 2*pi);
                    self.rbqz = mod(self.rbqz + randn(self.nRigidBodies, 1, 'single'), 2*pi);
                    self.rbqw = mod(self.rbqw + randn(self.nRigidBodies, 1, 'single'), 2*pi);
                    self.rbError = randn(self.nRigidBodies, 1, 'single');
                    self.rbTracked = uint8(rand(self.nRigidBodies, 1) > 0.02);
                    self.mx = self.mx + rand(self.nMarkers, 1, 'single');
                    self.my = self.my + rand(self.nMarkers, 1, 'single');
                    self.mz = self.mz + rand(self.nMarkers, 1, 'single');
                    self.msz = rand(self.nMarkers, 1, 'single');
                    self.mres = rand(self.nMarkers, 1, 'single') * 10e-4;
                else
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Rigid bodies

                    for r = 1:self.nRigidBodies
                        rb = frame.RigidBodies(r);
                        % If no rigid body exists in the current frame, rb may
                        % be an empty matrix. If so, create a struct with NaN
                        % fields
                        if isempty(rb)
                            rb = struct('x', nan, 'y', nan, 'z', nan, ...
                                'qx', nan, 'qy', nan, 'qz', nan, 'qw', nan, ...
                                'MeanError', nan, 'Tracked', 0);
                        end
                        self.rbx(r) = -rb.x;
                        self.rby(r) = rb.y;
                        self.rbz(r) = rb.z;
                        self.rbqx(r) = rb.qx;
                        self.rbqy(r) = rb.qy;
                        self.rbqz(r) = rb.qz;
                        self.rbqw(r) = rb.qw;
                        self.rbError(r) = rb.MeanError;
                        self.rbTracked(r) = uint8(rb.Tracked);

                        % Check if rigid body was successfully tracked. Issue
                        % warnings whenever it is lost, or on resuming after a lost
                        % period.
                        if self.rbTracked(r)
                            % Check if the last successfully tracked frame was
                            % the last frame. If not, this means we're resuming
                            % good tracking after a lost period. Notify the user
                            % with a warning
                            if self.lastTrackedFrameIdx(r) ~= self.frameIdx
                                eventData = RbEventData(r);
                                self.notify('trackingResumed', eventData);
                            end
                            self.lastTrackedFrameIdx(r) = newFrameIdx;
                        else
                            if self.lastTrackedFrameIdx(r) == self.frameIdx
                                eventData = RbEventData(r);
                                self.notify('trackingLost', eventData);
                            end
                            eventData = RbEventData(r);
                            self.notify('untrackedFrame', eventData);
                            self.nUntrackedFramesSinceWarn(r) = ...
                                self.nUntrackedFramesSinceWarn(r)+1;
                        end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Markers

                    for m = 1:self.nMarkers
                        marker = frame.LabeledMarkers(m);
                        if isempty(marker)
                            self.mx(m) = nan;
                            self.my(m) = nan;
                            self.mz(m) = nan;
                            self.msz(m) = nan;
                            self.mres(m) = nan;
                        else
                            self.mx(m) = marker.x;
                            self.my(m) = marker.y;
                            self.mz(m) = marker.z;
                            self.msz(m) = marker.size;
                            self.mres(m) = marker.residual;
                        end
                    end
                end
                
                self.nFramesAcquired = self.nFramesAcquired+1;
                self.nTotalFramesSinceWarn = self.nTotalFramesSinceWarn + 1;
                self.frameIdx = newFrameIdx;
                self.frameTimestamp = newFrameTime;
                if self.firstFrame, self.firstFrame = false; end
                
                self.logger.v('Frame #%u, Pos (m): [%.3f, %.3f, %.3f], Rot (rad): [%.3f %.3f %.3f %.3f]', ...
                    newFrameIdx, self.rbx, self.rby, self.rbz, self.rbqx, self.rbqy, self.rbqz, self.rbqw);
                
                if self.writeToFile
                    self.writeBuffer = [ ...
                        self.writeBuffer, self.currentFrameToBytes()];
                    
                    self.nFramesInBuffer = self.nFramesInBuffer + 1;
                    
                    if self.nFramesInBuffer >= self.writeBufferNFrames
                        self.flushBuffer();
                    end
                    
                end
                
                self.logger.v('Notifying "postGetFrame"');
                self.notify('postGetFrame');
                
            end
            
            % Update running average timer
            t = toc();
            self.meanGetFrameTime = runningMean( ...
                t, self.nFramesAcquired, self.meanGetFrameTime);
            self.logger.v('getFrame() execution time = %.1f ms\n', t*1e3);
            
        end
        
        function bytes = currentFrameToBytes(self)
            % CURRENTFRAMETOBYTES generate byte array from the current
            % frame of data
            bytes = [ ...
                typecast(self.frameIdx, 'uint8') ...                                    % int32 (4 bytes)
                typecast(self.frameTimestamp, 'uint8') ...                              % double (8 bytes)
                typecast([
                self.frameLatency
                self.rbx
                self.rby
                self.rbz
                self.rbqx
                self.rbqy 
                self.rbqz
                self.rbqw
                self.rbError]', 'uint8'), ...  % single (9*4 bytes)
                uint8(self.rbTracked)', ...% logical (1 byte)
                typecast([self.mx; self.my; self.mz; self.msz; self.mres], 'uint8')'
                ];
        end
        
        function printDescriptions(self, dataDescriptions)
            % Print out a description of actively tracked models from Motive
            
            lg = @(varargin) self.logger.i(varargin{:});
            
            lg('')
            lg('----------------------------------------------------------')
            lg('NatNet Tracking Models : %d', dataDescriptions.Count);
            
            for idx = 1 : dataDescriptions.Count
                descriptor = dataDescriptions.Item(idx-1);
                if(descriptor.type == 0)
                    lg('\tMarkerSet \t: %s', char(descriptor.Name));
                elseif(descriptor.type == 1)
                    lg('\tRigid Body \t: %s', char(descriptor.Name));
                elseif(descriptor.type == 2)
                    lg('\tSkeleton \t:  %n', char(descriptor.Name));
                else
                    lg('\tUnknown data type : ');
                end
            end
            
            lg('')
            
            for idx = 1 : dataDescriptions.Count
                descriptor = dataDescriptions.Item(idx-1);
                if(descriptor.type == 0)
                    lg('\tMarkerset : %s\t(%d markers)', ...
                        char(descriptor.Name), descriptor.nMarkers);
                    markerNames = descriptor.MarkerNames;
                    for markerIndex = 1 : descriptor.nMarkers
                        name = markerNames(markerIndex);
                        lg('\t\tMarker : %-20s\t(ID=%d)', char(name), ...
                            markerIndex);
                    end
                elseif(descriptor.type == 1)
                    lg('\tRigid Body : %s\t\t(ID=%d, ParentID=%d)', ...
                        char(descriptor.Name),descriptor.ID,descriptor.parentID);
                elseif(descriptor.type == 2)
                    lg('\tSkeleton : %s\t(%d bones)', char(descriptor.Name), descriptor.nRigidBodies);
                    %fprintf('\t\tID : %d\n', descriptor.ID);
                    rigidBodies = descriptor.RigidBodies;
                    for boneIndex = 1 : descriptor.nRigidBodies
                        rigidBody = rigidBodies(boneIndex);
                        lg('\tBone : %-20s\t(ID=%d, ParentID=%d)', ...
                            char(rigidBody.Name), rigidBody.ID, rigidBody.parentID);
                    end
                end
            end
            
            lg('');
            
        end
        
    end
    
    methods (Access = protected)
        function props = getPropertyGroups(self)
            import matlab.mixin.util.PropertyGroup
            props(1) = PropertyGroup({
                'hostIP'
                'fileName'
                'writeToFile'
                'frameIncrement'
                'nRigidBodies'
                'nMarkers'
                'simulate'
                'noDataTimeout'
                'untrackedWarnTreshold'
                'untrackedWarnInterval'}, 'User-defined settings');
            
            props(2) = PropertyGroup({
                'rbx'
                'rby'
                'rbz'
                'rbqx'
                'rbqy'
                'rbqz'
                'rbqw'
                'rbError'
                'rbTracked'
                'frameTimestamp'
                'frameIdx'}, 'Rigid body data');
            
            props(3) = PropertyGroup({
                'mx'
                'my'
                'mz'
                'msz'
                'mres'}, 'Marker data');
            
            props(4) = PropertyGroup({
                'frameRate'
                'nFramesAcquired'
                'nFramesDropped'
                'nFramesInBuffer'
                'meanGetFrameTime'
                'firstFrameIdx'
                'firstFrameTimestamp'}, 'Tracking summary data');
            
            props(5) = PropertyGroup({
                'fileOpen'
                'dataDir'
                'nWriteCycles'
                'meanWriteTime'
            }, 'Data file info');
        
        end
    end
    
end

function xb = runningMean(x, n, xb0)
% Incorporate new observation into a running mean
% x   - new observation
% n   - index of the new observation
% xb0 - last-calculated mean value

if n == 1
    xb = x;
else
    xb = (xb0*(n-1) + x) / n;
end

end
