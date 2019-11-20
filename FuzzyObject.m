classdef FuzzyObject < handle
    %FUZZYDEC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LocalSS
        MembershipFunctions
        ss_number
        MF_number
        x_ps
        x_history
        u_history
        u_ps
        k
    end
    
    methods
        function obj = FuzzyObject(LocalSS,MembershipFunctions, x_ps, u_ps)
            %FUZZYDEC Construct an instance of this class
            %   Detailed explanation goes here
            obj.LocalSS = LocalSS;
            obj.MembershipFunctions = MembershipFunctions;
            obj.ss_number = size(LocalSS);
            obj.MF_number = size(MembershipFunctions);
            obj.x_ps = x_ps;
            obj.u_ps = u_ps;

            assert(((obj.ss_number(1,1) == obj.MF_number(1,1)) & (obj.ss_number(1,2) == obj.MF_number(1,2))), "Models number and Membership functions must be the same");            
        end       
        
        function [] = reset(obj, x_0, u_0)
            obj.k = 0;
            obj.x_history = x_0;
            obj.u_history = u_0;
            
        end
        
        function [y,x] = countValue(obj,u)
            num = [0;0];
            den = [0;0];
            obj.k = obj.k+1;
            obj.u_history = [obj.u_history u];
            
            %Simulate all local objects.
            for i=1:obj.ss_number(1,1)
               for j= 1:obj.ss_number(1,2)
                   if obj.k > obj.LocalSS{i,j}.InputDelay(1,1)
                       v = obj.MembershipFunctions{i,j}.getValue(obj.u_history(1,obj.k - obj.LocalSS{i,j}.InputDelay(1,1)));
                       x_local = obj.LocalSS{i,j}.A*(obj.x_history(:,end)-obj.x_ps{i,j}) + obj.LocalSS{i,j}.B*(obj.u_history(:,obj.k - obj.LocalSS{i,j}.InputDelay(1,1))-obj.u_ps{i,j});
                       num  = num + v.*(x_local(:) + obj.x_ps{i,j});
                       den = den + v;
                   else
                       v = obj.MembershipFunctions{i,j}.getValue(obj.u_history(1));
                       num = num + v.*obj.x_ps{i,j};
                       den = den + v;
                   end
               end
            end
            
            %Global state of object.
            x = num./den;
            obj.x_history = [obj.x_history x];
            y = x(2,1);
            
        end
    end
end

