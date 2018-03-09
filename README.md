# matmot
This package provides some simple tools for streaming rigid body and marker data from Motive to disk, using OptiTrack's NatNet MATLAB client.

### Tasks handled
- Streaming from Motive to disk (Streamer.m)
- Loading previously streamed data into MATLAB (loadMtvFile.m)
- Identifying LEDs pulses for synchronization (extractSyncPulses.m)

### Requirements/dependencies
- MATLAB r2017a or later
- [MATLAB logging framework](https://github.com/brykko/matlab-logging)
- NaturalPoint Motive software (tested on v2, but should also work with v1)

### Setup
- Clone this repository and place it on your MATLAB path
- Install the MATLAB logging framework
- Configure your Motive software to broadcast tracking data in "multicast" mode

### Quick example script
Once you've followed the setup instructions above, you can acquire data as demonstrated by this example code snippet:
```matlab
% Create a Streamer object and acquire 60 seconds of data from 2 rigid bodies and 20 markers
streamer = matmot.Streamer('nRigidBodies', 2, 'nMarkers', 20);
streamer.start();
pause(60);
streamer.finish();
```
After finishing streaming, the data on disk can be loaded into MATLAB with the function `loadMtvFile()`.
