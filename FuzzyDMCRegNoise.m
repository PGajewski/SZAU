classdef FuzzyDMCRegNoise < handle
    %FUZZYDEC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LocalDMCs
        MembershipFunctions
        DMCs_number
        MF_number
        u_prev
        z_prev
    end
    
    methods
        function obj = FuzzyDMCRegNoise(LocalDMCs,MembershipFunctions)
            %FUZZYDEC Construct an instance of this class
            %   Detailed explanation goes here
            obj.LocalDMCs = LocalDMCs;
            obj.MembershipFunctions = MembershipFunctions;
            obj.DMCs_number = size(LocalDMCs);
            obj.MF_number = size(MembershipFunctions);
            assert(((obj.DMCs_number(1,1) == obj.MF_number(1,1)) & (obj.DMCs_number(1,2) == obj.MF_number(1,2))), "Models number and Membership functions must be the same");
        end
        
        function [] = reset(obj,u_p, Fdp)
            for i = 1:obj.DMCs_number(1,1)
                for j = 1:obj.DMCs_number(1,2)
                    obj.LocalDMCs{i,j}.reset(u_p, Fdp);
                    obj.u_prev = u_p;
                    obj.z_prev = Fdp;
                end
            end
        end
        
        function [] = setValue(obj, yzad)
            for i = 1:obj.DMCs_number(1,1)
                for j = 1:obj.DMCs_number(1,2)
                    obj.LocalDMCs{i,j}.setValue(yzad);
                end
            end
        end        
        
        function u = countValue(obj,y, z)
            num = 0;
            u_local = 0;
            for i=1:obj.DMCs_number(1,1)
               for j= 1:obj.DMCs_number(1,2)
                   v = obj.MembershipFunctions{i,j}.getValue(y);
                   u_local = u_local + obj.LocalDMCs{i,j}.countValue(y, z).*v;
                   num = num + v;
               end
            end
            %u = num/den;
            u = u_local./num;
            %Update regulator for real u value.
            for i=1:obj.DMCs_number(1,1)
               for j= 1:obj.DMCs_number(1,2)
                   obj.LocalDMCs{i,j}.u_prev=u;
                   obj.LocalDMCs{i,j}.deltaup(1) = u - obj.u_prev;
                   obj.LocalDMCs{i,j}.z_prev=z;
                   obj.LocalDMCs{i,j}.deltazp(1) = z - obj.z_prev;
               end
            end
            obj.u_prev = u;
            obj.z_prev = z;
        end
    end
end

