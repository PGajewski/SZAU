function LinearPlotsEx1(F1JumpVal, F1JumpTime, FdJumpVal, FdJumpTime, sim_time)
% LinearPlotsEx1(10, 500, 10, 200, 2000)
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

F1=[0, 0;
    F1JumpTime, F1JumpVal];
FD=[0, 0;
    FdJumpTime, FdJumpVal];


[t,h] = linearSimulation(F1,FD,sys, Ts, sim_time, [0,0]);

h = h + ones(size(h)).*[h1p,h2p];

F1vector = ones(sim_time, 1);
F1vector(1:F1(2,1)) = F1(1,2)+F1p;
F1vector(F1(2,1):sim_time) = F1(2,2)+F1p;

Fdvector = ones(sim_time, 1);
Fdvector(1:FD(2,1)) = FD(1,2)+FDp;
Fdvector(FD(2,1):sim_time) = FD(2,2)+FDp;

subplot(2,1,1);
plot(t,h);
legend('h1','h2 - wyj�cie', 'Location', 'northwest');
xlabel('Czas [s]');
ylabel('Stany wewnetrzne');
title('Przebieg zmiennych stanu po linearyzacji');

subplot(2,1,2);
plot(F1vector);
hold on;
plot(Fdvector);
hold off;
legend('F1','Fd - zak��cenie', 'Location', 'west');
xlabel('Czas [s]');
ylabel('Wej�cia obiektu');
title(sprintf('Przebieg wej�� po linearyzacji'));
%print(sprintf('pdfs/SymulacjaLiniowaF1=%dFd=%d', F1JumpVal, FdJumpVal), '-dpdf');
