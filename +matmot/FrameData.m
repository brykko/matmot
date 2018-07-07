classdef FrameData < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frameIdx
        frameTimestamp
        frameLatency
        rbx
        rby
        rbz
        rbqx
        rbqy
        rbqz
        rbqw
        rbError
        rbTracked
        mx
        my
        mz
        msz
        mres
    end
    
    methods
        function bytes = toBytes(self)
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
        
        
    end
    
end

