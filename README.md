# matmot
This package provides some simple tools for streaming rigid body and marker data from Motive to disk, using OptiTrack's NatNet MATLAB client.

### Tasks handled
- Streaming from Motive to disk (MotiveStreamer.m)
- Loading previously streamed data into MATLAB (loadMtvFile.m)
- Identifying LEDs pulses for synchronization (extractMotiveSyncPulses.m)

### Requirements/dependencies
- Matlab r2017a or later
- [Matlab logging framework](https://github.com/brykko/matlab-logging)
