classdef RbEventData < event.EventData
    %RBEVENTDATA conveys information about a rigid body event.
    %
    % This subclass of EventData adds the "rigidBodyId" property which
    % allows event listeners to access the identity of the rigid body the
    % event refers to. The streaming ID of a rigid body is configured in
    % Motive.
    
    properties
        rigidBodyId
    end
    
    methods
        function self = RbEventData(rigidBodyId)
            self.rigidBodyId = rigidBodyId;
        end
    end
    
end

