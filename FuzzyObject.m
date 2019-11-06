classdef FuzzyObject < handle
    %FUZZYDEC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LocalTfs
        MembershipFunctions
        Tf_number
        MF_number
        x_ps
        x_prev
        u_ps
    end
    
    methods
        function obj = FuzzyObject(LocalTfs,MembershipFunctions, x_ps, u_ps)
            %FUZZYDEC Construct an instance of this class
            %   Detailed explanation goes here
            obj.LocalTfs = LocalTfs;
            obj.MembershipFunctions = MembershipFunctions;
            obj.Tf_number = size(LocalTfs);
            obj.MF_number = size(MembershipFunctions);
            obj.x_ps = x_ps;
            obj.u_ps = u_ps;
            
            % Clear all delays.
            for i=1:obj.Tf_number(1,1)
               for j= 1:obj.Tf_number(1,2)
                   s = size(obj.LocalTfs{i,j}.InputDelay)
                   obj.LocalTfs{i,j}.InputDelay = zeros(s);
               end
            end
            
            assert(((obj.Tf_number(1,1) == obj.MF_number(1,1)) & (obj.Tf_number(1,2) == obj.MF_number(1,2))), "Models number and Membership functions must be the same");            
        end       
        
        function [] = reset(obj, x_0)
            obj.x_prev = x_0;
            
        end
        
        function [y,x] = countValue(obj,u)
            num = [0 0];
            den = [0 0];
            %Simulate all local objects.
            for i=1:obj.Tf_number(1,1)
               for j= 1:obj.Tf_number(1,2)
                   v = obj.MembershipFunctions{i,j}.getValue(u(1));
                   [y_local, t, x_local] = lsim(obj.LocalTfs{i,j},[u-obj.u_ps{i,j} u-obj.u_ps{i,j}],[],obj.x_prev-obj.x_ps{i,j});
                   num  = num + v.*(x_local(end,:) + obj.x_ps{i,j});
                   den = den + v;
               end
            end
            
            %Global state of object.
            x = num./den;
            obj.x_prev = x;
            y = x(1,2);
            
        end
    end
end

