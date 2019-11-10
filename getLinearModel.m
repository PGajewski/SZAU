function [h] = getLinearModel(F1p, FDp)
%GETLINEARMODEL Summary of this function goes here
%   Detailed explanation goes here
% FDp is 14 by default

alfa1 = 20;
alfa2 = 22;

h = zeros(2,1);
h(1) = ((F1p+FDp)/alfa1)^2;
h(2) = ((F1p+FDp)/alfa2)^2;
end

