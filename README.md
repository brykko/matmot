# matmot
This package provides some simple tools for streaming rigid body and marker data from Motive to disk, using OptiTrack's NatNet MATLAB client.

### Tasks handled
- Streaming from Motive to disk (Streamer.m)
- Loading previously streamed data into MATLAB (loadMtvFile.m)
- Identifying LEDs pulses for synchronization (extractSyncPulses.m)

### Requirements/dependencies
- MATLAB r2017a or later
- [MATLAB logging framework](https://github.com/brykko/matlab-logging)

### Setup
- Clone this repository and place it on your MATLAB path
- Install the MATLAB logging framework
