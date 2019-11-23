function [t,x,y] = objectSimulation(F1_with_time,FD_with_time, tau, sim_time, x_0)
%objectSimulation Simulation of object

% Check variables dimension.
temp1 = size(F1_with_time);
temp2 = size(FD_with_time);
assert((temp1(2) == 2 && temp2(2) == 2),'Incompatible changing moment values!');
assert((isnumeric(tau)),'Transport delay must be a number!');
assert((isvector(x_0) && length(x_0) == 2),'Init condition must be a vector with size 2'); 

%Check last value change moment timestep.
assert((F1_with_time(end,1)<=sim_time) && (FD_with_time(end,1)<=sim_time),'One or more value changes is outside simulation range!');

%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%Add simulation end to changes moments.
F1 = [F1_with_time;[sim_time,F1_with_time(end,2)]];
FD = [FD_with_time;[sim_time,FD_with_time(end,2)]];

% Merge all changing signal variable into one (ugly version).
% Init.
F1_length = temp1(1)+1;
FD_length = temp2(1)+1;
counter_F1 = 2;
counter_FD = 2;


values_with_time = [0;F1(1,2);FD(1,2)];
actual_time = min([F1(2,1) FD(2,1)]);

while 1
    %Add next timestamp with values.
    actual_F1 = 0;
    actual_FD = 0;
    %Increase counters.
    if (actual_time ~= tau) && (actual_time >= (F1(counter_F1,1)+tau)) && (counter_F1 ~= F1_length)
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

    values_with_time = [values_with_time [actual_time;actual_F1;actual_FD]];
    
    if actual_time == sim_time
        break;
    end
     %Find minimal time.
    actual_time = min([(F1(counter_F1,1) + tau) (FD(counter_FD,1))]);
end

% Main simulation.
x = [];
t = [];

s = size(values_with_time);
init_values = x_0;
% Simulate by all changes moments.
% u = [];
for i=1:(s(2)-1)
    stateHandler = @(t,x) stateFunction(t,x,values_with_time(2,i), values_with_time(3,i));
    [temp1,temp2]=ode45(stateHandler,[values_with_time(1,i) values_with_time(1,i+1)],init_values, options);
    %s_t = size(temp1);
    t = [t temp1(1:end-1,:)'];
    x = [x temp2(1:end-1,:)'];
    %u = [u [ones(1,s_t(1)-1)*values_with_time(2,i);ones(1,s_t(1)-1)*values_with_time(3,i)]];
    % Remember previous state for next step.
    init_values = temp2(end,:);
end

t =[t sim_time];
x =[x init_values'];
%u = [u [values_with_time(2,end);values_with_time(3,end)]];
process_length = length(t);

%Count output.
%Hight.
y(1,:) = x(2,:);
end

