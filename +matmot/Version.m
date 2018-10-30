classdef Version < matlab.mixin.CustomDisplay
    %VERSION simple class for comparing sequential versions (e.g. x.x.x)
    
    properties
        version
    end
    
    methods
        function self = Version(ver)
            %VERSION
            if ischar(ver)
                self = matmot.Version.fromString(ver);
            elseif isnumeric(ver)
                self.version = ver;
            elseif isa(ver, 'matmot.Version')
                self.version = ver.version;
            end
        end
        
        function str = toString(self)
            v = self.version;
            str = strjoin(strsplit(num2str(v), ' '), '.');
        end
        
        function verDiff = difference(self, ver)
            verObj = matmot.Version(ver);
            v1 = self.version;
            v2 = verObj.version;
            assert(numel(v1)==numel(v2), 'matmot:Version:unequalLength', ...
                'Versions must have equal number of elements to compare.')
            verDiff = v1-v2;
        end
        
        function tf = gt(self, ver)
            verDiff = self.difference(ver);
            idx = find(verDiff ~= 0, 1);
            tf = ~isempty(idx) && verDiff(idx) > 0;
        end
        
        function tf = ge(self, ver)
            tf = self.gt(ver) || self.eq(ver);
        end
        
        function tf = lt(self, ver)
            verDiff = self.difference(ver);
            idx = find(verDiff ~= 0, 1);
            tf = ~isempty(idx) && verDiff(idx) < 0;
        end
        
        function tf = le(self, ver)
            tf = self.lt(ver) || self.eq(ver);
        end
        
        function tf = eq(self, ver)
            verDiff = self.difference(ver);
            idx = find(verDiff ~= 0, 1);
            tf = isempty(idx);
        end
    end
    
    methods (Access=protected)
        function displayScalarObject(self)
            fprintf('%s\n\n', self.toString());
        end
    end
    
    methods (Static)
        function obj = fromString(verStr)
            verElStr = strsplit(verStr, '.');
            for i = 1:numel(verElStr)
                ver(i) = str2num(verElStr{i});
            end
            obj = matmot.Version(ver);
        end
    end
end
