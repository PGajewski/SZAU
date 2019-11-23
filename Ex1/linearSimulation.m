function [t,x,y] = linearSimulation(F1_with_time,FD_with_time, G, Ts, sim_time, x_0)
%objectSimulation Simulation of object
%   FH_with_time - matrix of changing FH values: new value (first column) and time step (second column).
%   FC_with_time - matrix of changing FC values: new value (first column) and time step (second column).
%   FD_with_time - matrix of changing FD values: new value (first column) and time step (second column).
%   tau_H        - transport delay of signal FH.
%   tau_C        - transport delay of signal FC.
%   tau_D        - transport delay of signal FD.

% Check variables dimension.
temp1 = size(F1_with_time);
temp2 = size(FD_with_time);
assert((temp1(2) == 2 && temp2(2) == 2),'Incompatible changing moment values!');
assert((isvector(x_0) && length(x_0) == 2),'Init condition must be a vector with size 2'); 

%Check last value change moment timestep.
assert((F1_with_time(end,1)<=sim_time) && (FD_with_time(end,1)<=sim_time),'One or more value changes is outside simulation range!');

F1_length = temp1(1)+1;
FD_length = temp2(1)+1;

counter_F1 = 2;
counter_FD = 2;

%Add simulation end to changes moments.
F1 = [F1_with_time;[sim_time,F1_with_time(end,2)]];
FD = [FD_with_time;[sim_time,FD_with_time(end,2)]];

t = 0:Ts:sim_time;
values = zeros(length(t),2);
values(1,:) = [F1_with_time(1,2),FD_with_time(1,2)];
% Merge all changing signal variable into one (ugly version).
for i=1:length(t)
% Main simulation.
    %Increase counters.
    actual_time = t(i);
    actual_F1 = 0;
    actual_FC = 0;

    %Increase counters.
    if (actual_time >= F1(counter_F1,1)) && (counter_F1 ~= F1_length)
        actual_F1 = F1(counter_F1,2);
        counter_F1 = counter_F1 + 1;
    else
        actual_F1 = F1(counter_F1-1,2);
    end

    if (actual_time >= FD(counter_FD,1)) && (counter_FD ~= FD_length)
        actual_FD = FD(counter_FD,2);
        counter_FD = counter_FD + 1;
    else
        actual_FD = FD(counter_FD-1,2);
    end

    values(i,:) = [actual_F1 actual_FD];
    
    if actual_time == sim_time
        break;
    end
     %Find minimal time.
end
[y,t,x]=lsim(G, values, t ,x_0);
end

