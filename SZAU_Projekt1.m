clear all;
close all;

%wpisanie danych
F1p=73;
FDp=14;
tau=150;
h2p=15.6384;
h1p=18.9225;
Ts = 0.1;
C1 = 0.35;
C2 = 0.3;
alfa1 = 20;
alfa2 = 22;

%% Test work point.
F1=[0, F1p;
    500, F1p+1];
FD=[0, FDp;
    1200, FDp+1];

sim_time = 2000;
[t1,x1,y1] = objectSimulation(F1,FD,tau, sim_time, [h1p, h2p]);

figure(1);
plot(t1,x1);
legend('h1','h2');
xlabel('Czas [s]');
ylabel('Stany wewnetrzne');
title('Przebieg zmiennych stanu (symulacja)');

figure(2);
plot(t1,y1);
legend('h2')
xlabel('czas [s]');
ylabel('Wyjscia obiektu');
title('Przebieg wyjscia obiektu (symulacja)');

%% Linearization.
a11 = alfa1 * h1p^(-5/2)/(2*C1) - 2/(3*C1) * (F1p + FDp) * h1p^(-3);
a12 = 0;
a21 = alfa1/(6*C2) * h1p^(-1/2) * h2p^-2;
a22 = -2 * alfa1/(3*C2) * sqrt(h1p) * h2p^-3 + alfa2/(2*C2) * h2p^(-5/2);
u11 = h1p^(-2)/(3*C1);
u12 = h1p^(-2)/(3*C1);
u21 = 0;
u22 = 0;

A = [a11, a12; a21, a22];
B = [u11 u12; u21 u22];
C = [0 1];
D = [0 0];

sys = ss(A,B,C,D,'InputDelay',[tau,0]);
%sys = ss(A,B,C,D);

%% Linear simulation
F1=[0, 0;
    500, 1];
FD=[0, 0;
    1200, 1];

[t2,x2,y2] = linearSimulation(F1,FD, sys, Ts, sim_time, [0, 0]);

x2 = x2 + ones(size(x2)).*[h1p,h2p];
y2 = y2 + ones(size(y2)).*h2p;

figure(3);
plot(t2,x2);
legend('h1','h2');
xlabel('Czas [s]');
ylabel('Stany wewnetrzne');
title('Przebieg zmiennych stanu (obiekt zlinearyzowany)');

figure(4);
plot(t2,y2);
legend('h2')
xlabel('czas [s]');
ylabel('Wyjscia obiektu');
title('Przebieg wyjscia obiektu (obiekt zlinearyzowany)');

%% Compare responses od non-linear and linear model for different steps.
steps_number = 5;
step_value = 2;
% F1 step
figure(5);
for i=1:steps_number
    F1=[0, F1p;
        500, F1p+i*step_value];
    FD=[0, FDp];
    [t,x,y] = objectSimulation(F1,FD,tau, sim_time, [h1p, h2p]);
    plot(t,y);
    hold on;
    
    F1=[0, F1p;
        500, F1p-i*step_value];
    FD=[0, FDp];
    [t,x,y] = objectSimulation(F1,FD,tau, sim_time, [h1p, h2p]);
    plot(t,y);
    hold on;
end
for i=1:steps_number
    F1=[0, 0;
        500, i*step_value];
    FD=[0, 0];
    [t,x,y] = linearSimulation(F1,FD, sys, Ts, sim_time, [0, 0]);
    x = x + ones(size(x)).*[h1p,h2p];
    y = y + ones(size(y)).*h2p;
    plot(t,y);
    hold on;
    
    F1=[0, 0;
        500, -i*step_value];
    FD=[0, 0];
    [t,x,y] = linearSimulation(F1,FD, sys, Ts, sim_time, [0, 0]);
    x = x + ones(size(x)).*[h1p,h2p];
    y = y + ones(size(y)).*h2p;
    plot(t,y);
    hold on;
end

% FD step
figure(6);
for i=1:steps_number
    F1=[0, F1p];
    FD=[0, FDp;
        500, FDp+i*step_value];
    [t,x,y] = objectSimulation(F1,FD,tau, sim_time, [h1p, h2p]);
    plot(t,y);
    hold on;
    
    F1=[0, F1p];
    FD=[0, FDp;
        500, FDp-i*step_value];
    [t,x,y] = objectSimulation(F1,FD,tau, sim_time, [h1p, h2p]);
    plot(t,y);
    hold on;
end
for i=1:steps_number
    F1=[0, 0];
    FD=[0, 0;
        500, i*step_value];
    [t,x,y] = linearSimulation(F1,FD, sys, Ts, sim_time, [0, 0]);
    x = x + ones(size(x)).*[h1p,h2p];
    y = y + ones(size(y)).*h2p;    
    plot(t,y);
    hold on;
    
    F1=[0, 0];
    FD=[0, 0;
        500, -i*step_value];
    [t,x,y] = linearSimulation(F1,FD, sys, Ts, sim_time, [0, 0]);
    x = x + ones(size(x)).*[h1p,h2p];
    y = y + ones(size(y)).*h2p;
    plot(t,y);
    hold on;
end

%% Prediction DMC
T=0.5;
Tk = 4000;
%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%%Prepare discrete transfer function.
Gz = c2d(tf(ss(A,B,C,D)),T,'zoh');
Gz.InputDelay = [tau/T,0];

%%Prepare DMC regulator.
dmc = DMCReg(Gz,2393, 600, 1, 1);
dmc.reset(F1p);
dmc.setValue(h2p+0.5);

uk= ones((Gz.InputDelay(1)),1).*F1p;
y = ones(Tk, 1).*h2p;
h = [h1p, h2p];

%Main simulation loop.
for k=2:Tk
    if k > (Gz.InputDelay(1))
        stateHandler = @(t,x) stateFunction(t,x,uk(k - (Gz.InputDelay(1))), FDp);
        [t, h] = ode45(stateHandler,[0 Gz.Ts],h(end, :), options);
        y(k) = h(end,2);
    end
    uk(k) = dmc.countValue(y(k));
end

figure(7);
subplot(2,1,2);
stairs( uk, 'g');
%title(strcat('Dzialanie regulatora dla nastaw D=', num2str(D), ' N =', num2str(N), ' Nu=', num2str(Nu), ' lambda=',num2str(lambda)));
legend('Sterowanie');
subplot(2,1,1);
stairs(ones(Tk,1).*(h2p+0.5), 'b');
hold on;
stairs(y, 'r');
legend('Wyjœcie zadane', 'wyjœcie regulatora', 'Location', 'east');
xlabel('k');
ylabel('y/u');
hold off;

%% Different models for fuzzy logics.
F1p_array = [(F1p-60) (F1p-30) F1p (F1p+30) (F1p+60)];
h1p_array = [];
h2p_array = [];
sim_time = 2000;
tfs = cell(size(F1p_array));
sss = cell(size(F1p_array));
T=0.5;

i = 1;
for local_F1p=F1p_array
    %Simulate non linear object to get h vector of new work point.
    [t,x,y] = objectSimulation([0 local_F1p],[0 FDp],tau, sim_time, [0.1, 0.1]);
    x_p = x(:,end);
    h1p_array = [h1p_array x_p(1)];
    h2p_array = [h2p_array x_p(2)];
    %Linearization.
    a11 = alfa1 * x_p(1)^(-5/2)/(2*C1) - 2/(3*C1) * (local_F1p + FDp) * x_p(1)^(-3);
    a12 = 0;
    a21 = alfa1/(6*C2) * x_p(1)^(-1/2) * x_p(2)^-2;
    a22 = -2 * alfa1/(3*C2) * sqrt(x_p(1)) * x_p(2)^-3 + alfa2/(2*C2) * x_p(2)^(-5/2);
    u11 = x_p(1)^(-2)/(3*C1);
    u12 = x_p(1)^(-2)/(3*C1);
    u21 = 0;
    u22 = 0;

    A = [a11, a12; a21, a22];
    B = [u11 u12; u21 u22];
    C = [0 1];
    D = [0 0];

    tfs{i} = c2d(tf(ss(A,B,C,D,'InputDelay',[tau,0])),T,'zoh');
    tfs{i}.InputDelay = [tau/T,0];
    sss{i} = c2d(ss(A,B,C,D),T,'zoh');

    i = i + 1;
end

%% Fuzzy object.
T=0.5;
Tk = 1500;

%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%Prepare local DMCs.
LocalObjects = cell(1,5);
LocalObjects{1}=sss{1};
LocalObjects{2}=sss{2};
LocalObjects{3}=sss{3};
LocalObjects{4}=sss{4};
LocalObjects{5}=sss{5};

%Prepare membership functions.
MembershipFunctions = cell(1,5);
MembershipFunctions{1}=MembershipFunction([(F1p-50) 1; (F1p-40) 0]);
MembershipFunctions{2}=MembershipFunction([(F1p-50) 0; (F1p-40) 1; (F1p-20) 1; (F1p-10) 0]);
MembershipFunctions{3}=MembershipFunction([(F1p-20) 0; (F1p-10) 1; (F1p+10) 1; (F1p+20) 0]);
MembershipFunctions{4}=MembershipFunction([(F1p+10) 0; (F1p+20) 1; (F1p+40) 1; (F1p+50) 0]);
MembershipFunctions{5}=MembershipFunction([(F1p+40) 0; (F1p+50) 1]);

%Prepare x_0 of models.
x_0_models = cell(1,5);
x_0_models{1} = [h1p_array(1) h2p_array(1)];
x_0_models{2} = [h1p_array(2) h2p_array(2)];
x_0_models{3} = [h1p_array(3) h2p_array(3)];
x_0_models{4} = [h1p_array(4) h2p_array(4)];
x_0_models{5} = [h1p_array(5) h2p_array(5)];

%Prepare x_0 of models.
u_0_models = cell(1,5);
u_0_models{1} = [F1p_array(1); FDp];
u_0_models{2} = [F1p_array(2); FDp];
u_0_models{3} = [F1p_array(3); FDp];
u_0_models{4} = [F1p_array(4); FDp]; 
u_0_models{5} = [F1p_array(5); FDp];

%Main simulation.
%%Prepare DMC regulator.
fuzzy_obj = FuzzyObject(LocalObjects, MembershipFunctions,x_0_models, u_0_models);
fuzzy_obj.reset([h1p h2p]);

uk= ones(Tk,1).*(F1p+5);
y1 = ones(Tk, 1).*h2p;
y2 = ones(Tk, 1).*h2p;

h1 = [h1p, h2p];
h2 = [h1p, h2p];

%Main simulation loop.
for k=2:Tk
    if k > tau/T
        %Non linear object.
        stateHandler = @(t,x) stateFunction(t,x,uk(k - (tau/T)), FDp);
        [t, h] = ode45(stateHandler,[0 T],h1(end, :), options);
        h1 = [h1; h];
        y1(k) = h1(end,2);
        
        %Fuzzy object.
        [temp1, temp2] = fuzzy_obj.countValue([uk(k - tau/T);FDp]);
        y2(k) = temp1;
        h2 = [h2; temp2];
    end
    k
end
figure(8);
plot(T:T:Tk*T,y1);
hold on;
plot(T:T:Tk*T,y2);
title('Porownanie obiektu nieliniowego oraz rozmytego');

%% Fuzzy DMC.
T=0.5;
Tk = 5000;
%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%Prepare local DMCs.
LocalDMCs = cell(1,5);
LocalDMCs{1}=DMCReg(tfs{1},2393, 600, 1, 1);
LocalDMCs{2}=DMCReg(tfs{2},2393, 600, 1, 1);
LocalDMCs{3}=DMCReg(tfs{3},2393, 600, 1, 1);
LocalDMCs{4}=DMCReg(tfs{4},2393, 600, 1, 1);
LocalDMCs{5}=DMCReg(tfs{5},2393, 600, 1, 1);

%Prepare membership functions.
MembershipFunctions = cell(1,5);
MembershipFunctions{1}=MembershipFunction([(h2p_array(1)+2) 1; (h2p_array(2)-2) 0]);
MembershipFunctions{2}=MembershipFunction([(h2p_array(1)+2) 0; (h2p_array(2)-2) 1; (h2p_array(2)+2) 1; (h2p_array(3)-2) 0]);
MembershipFunctions{3}=MembershipFunction([(h2p_array(2)+2) 0; (h2p_array(3)-2) 1; (h2p_array(3)+2) 1; (h2p_array(4)-2) 0]);
MembershipFunctions{4}=MembershipFunction([(h2p_array(3)+2) 0; (h2p_array(4)-2) 1; (h2p_array(4)+2) 1; (h2p_array(5)-2) 0]);
MembershipFunctions{5}=MembershipFunction([(h2p_array(4)+2) 0; (h2p_array(5)-2) 1]);

%Main simulation.
%%Prepare DMC regulator.
fuzzy_dmc = FuzzyDMCReg(LocalDMCs, MembershipFunctions);
fuzzy_dmc.reset(F1p);
fuzzy_dmc.setValue(h2p+5);

uk= ones((LocalDMCs{1}.Gz.InputDelay(1)),1).*F1p;
y = ones(Tk, 1).*h2p;
h = [h1p, h2p];

%Main simulation loop.
for k=2:Tk
    if k > (LocalDMCs{1}.Gz.InputDelay(1))
        stateHandler = @(t,x) stateFunction(t,x,uk(k - (LocalDMCs{1}.Gz.InputDelay(1))), FDp);
        [t, h] = ode45(stateHandler,[0 LocalDMCs{1}.Gz.Ts],h(end, :), options);
        y(k) = h(end,2);
    end
    uk(k) = fuzzy_dmc.countValue(y(k));
    k
end

figure(9);
subplot(2,1,2);
stairs( uk, 'g');
%title(strcat('Dzialanie regulatora dla nastaw D=', num2str(D), ' N =', num2str(N), ' Nu=', num2str(Nu), ' lambda=',num2str(lambda)));
legend('Sterowanie');
subplot(2,1,1);
stairs(ones(Tk,1).*(h2p+5), 'b');
hold on;
stairs(y, 'r');
legend('Wyjœcie zadane', 'wyjœcie regulatora', 'Location', 'east');
xlabel('k');
ylabel('y/u');
hold off;

