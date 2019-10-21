clear all;
close all;

%wpisanie danych
C=0.4;
a=8;
F1p=73;
FDp=14;
tau=150;
h1p=15.6384;
h2p=18.9225;
Ts = 0.1;


%% Test work point.
F1=[0, F1p];
FD=[0, FDp];

sim_time = 2000;
[t1,x1,y1] = objectSimulation(F1,FD,tau, sim_time, [0.01, 0.01]);

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

h1p = x1(1,end);
h2p = x1(2,end);

%%Linearization.

