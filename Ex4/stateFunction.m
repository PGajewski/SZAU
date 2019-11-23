function [h_next] = stateFunction(t,h,F_1, F_D)
%stateFunction - equation of state in model.
C1 = 0.35;
C2 = 0.3;
alfa1 = 20;
alfa2 = 22;
h_1 = (1/(3*C1*h(1)^2))*(F_1 + F_D - alfa1*sqrt(h(1)));
h_2 = (1/(3*C2*h(2)^2))*(alfa1*sqrt(h(1)) - alfa2*sqrt(h(2)));
h_next = [h_1;h_2];
end

