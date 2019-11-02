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
