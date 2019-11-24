classdef FuzzyObjectK < handle
    %FUZZYDEC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        LocalModels
        MembershipFunctions
        ss_number
        MF_number
    end
    
    methods
        function obj = FuzzyObjectK(LocalModels,MembershipFunctions)
            obj.LocalModels = LocalModels;
            obj.MembershipFunctions = MembershipFunctions;
            obj.ss_number = size(LocalModels);
            obj.MF_number = size(MembershipFunctions);
        end       
        
        function [y, x] = countValue(obj,u)
            x = [0,0];
            C1 = 0.35;
            C2 = 0.3;
            alfa1 = 20;
            alfa2 = 22;
            Fdp = 14;
            sumMF = 0; % norm
            %Simulate all local objects.
            for i=1:obj.ss_number(1,1)
               for j= 1:obj.ss_number(1,2)
                   hp = getLinearModel(obj.LocalModels(j), 14);
                   h1 = @(x) (-5*alfa1/(6*C1)*hp(1)^(-3/2)+2*hp(1)^(-2)*(obj.LocalModels(j)+Fdp)/(3*C1)+hp(1)^(-2)*(x+Fdp)/(3*C1))/-(alfa1*hp(1)^(-5/2)/(2*C1)-2*(obj.LocalModels(j)+Fdp)*hp(1)^(-3)/(3*C1));
                   h2 = @(x) (5*alfa1*hp(1)^(1/2)*hp(2)^(-2)/(6*C2)-5*alfa2*hp(2)^(-3/2)/(6*C2)+alfa1*hp(1)^(-1/2)*hp(2)^(-2)*h1(x)/(6*C2))/(2*alfa1*hp(1)^(1/2)*hp(2)^(-3)/(3*C2)-alfa2*hp(2)^(-5/2)/(2*C2));
                   x = x + [h1(u), h2(u)]*obj.MembershipFunctions{i,j}.getValue(u);
                   sumMF = sumMF + obj.MembershipFunctions{i,j}.getValue(u);
               end
            end
            x = x./sumMF;
            y=x(1,2);
        end
    end
end

