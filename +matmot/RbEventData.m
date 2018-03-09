classdef RbEventData < event.EventData
    %RBEVENTDATA conveys information about a rigid body event.
    %
    % This subclass of EventData adds the "rigidBodyId" property which
    % allows event listeners to access the identity of the rigid body the
    % event refers to. The ID is the index of the rigid body (1-based) as
    % broadcast by Motive and received by the NatNet client.
    
    properties
        rigidBodyId
    end
    
    methods
        function self = RbEventData(rigidBodyId)
            self.rigidBodyId = rigidBodyId;
        end
    end
    
end

