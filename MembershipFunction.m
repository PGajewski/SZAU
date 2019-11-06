classdef MembershipFunction < handle
    %MembershipFunction Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        s
        charPoints
    end
    
    methods
        function obj = MembershipFunction(charPoints)
            obj.s = size(charPoints);
            for i=2:obj.s(1,1)
                if charPoints(i,1)<=charPoints(i-1,1)
                   error("Characteristic points must be rising!");
                end
            end
            obj.charPoints = charPoints;
        end
        
        function output = getValue(obj,value)
            %characteristic points.

            %Border conditions.
            if obj.charPoints(1,1) >= value
                output = obj.charPoints(1,2);
            elseif obj.charPoints(end,1) <= value
                output = obj.charPoints(end,2);
            else
                %Find value
                for i=2:obj.s(1,1)
                    if obj.charPoints(i,1) > value
                        output = (obj.charPoints(i,2)-obj.charPoints(i-1,2))/(obj.charPoints(i,1)-obj.charPoints(i-1,1))*(value-obj.charPoints(i-1,1))+obj.charPoints(i-1,2);
                        break
                    elseif obj.charPoints(i,1) == value
                        output = obj.charPoints(i,2);
                    end
                end
            end
        end
    end
end

