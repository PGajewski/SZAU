classdef FuzzyDEC
    %FUZZYDEC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LocalDECs
        MembershipFunctions
    end
    
    methods
        function obj = FuzzyDEC(LocalDECs,MembershipFunctions)
            %FUZZYDEC Construct an instance of this class
            %   Detailed explanation goes here
            obj.LocalDECs = LocalDECs;
            obj.MembershipFunctions = MembershipFunctions;
        end
        
        function u = count(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            u = 0;
        end
    end
end

