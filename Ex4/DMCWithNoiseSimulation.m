load Gz.mat;
D = 2393;
N = 600;
Nu = 1;
lambda =1;
yzad = h2p+0.5;
Tk = 4000;
Fd = [ ones(1,Tk/2).*14, ones(1,Tk/2).*16 ] ;Tk = 4000;

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%%Prepare DMC regulator.
dmc = DMCRegWithNoise(Gz,D, N, Nu, lambda);
dmc.reset(F1p, FDp);
dmc.setValue(h2p+0.5);

uk= ones((Gz.InputDelay(1)),1).*F1p;
y = ones(Tk, 1).*h2p;
h = [h1p, h2p];

%Main simulation loop.
for k=2:Tk
    if k > (Gz.InputDelay(1))
        stateHandler = @(t,x) stateFunction(t,x,uk(k - (Gz.InputDelay(1))), Fd(k));
        [t, h] = ode45(stateHandler,[0 Gz.Ts],h(end, :), options);
        y(k) = h(end,2);
    end
    uk(k) = dmc.countValue(y(k), Fd(k));
end

figure();
subplot(2,1,1);
plot(ones(Tk,1).*(yzad), 'b-');
title(strcat('Regulator z uwzglêdnieniem zak³óceñ dla nastaw D=', num2str(D), ' N =', num2str(N), ' Nu=', num2str(Nu), ' lambda=',num2str(lambda)));
hold on;
stairs(y, 'r');
legend('Wyjœcie zadane', 'Wyjœcie regulatora', 'Location', 'east');
xlabel('k');
ylabel('y');
subplot(2,1,2);
stairs( uk, 'g');
legend('Sterowanie');
xlabel('k');
ylabel('u');
title(strcat('Skok wartoœci zadanej =', num2str(yzad-h2p), ' Skok zak³ócenia o ', num2str(Fd(end)-Fd(1))));
hold off;
%print(sprintf('DMCNoiseskok=%.1f.pdf', round(yzad-h2p, 1)), '-dpdf');