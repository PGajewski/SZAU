%% Prediction DMC with noise
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

T=0.5;
Tk = 4000;
%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%%Prepare discrete transfer function.
Gz = c2d(tf(ss(A,B,C,D)),T,'zoh');
Gz.InputDelay = [tau/T,0];

%%Prepare DMC regulator.
dmc = DMCRegWithNoise(Gz,2393, 600, 1, 1);
dmc.reset(F1p, FDp);
dmc.setValue(h2p+0.5);

uk= ones((Gz.InputDelay(1)),1).*F1p;
y = ones(Tk, 1).*h2p;
h = [h1p, h2p];
Fd = FDp+3/Tk:3/Tk:FDp+3;
%Main simulation loop.
for k=2:Tk
    if k > (Gz.InputDelay(1))
        stateHandler = @(t,x) stateFunction(t,x,uk(k - (Gz.InputDelay(1))), Fd(k));
        [t, h] = ode45(stateHandler,[0 Gz.Ts],h(end, :), options);
        y(k) = h(end,2);
    end
    uk(k) = dmc.countValue(y(k), Fd(k));
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