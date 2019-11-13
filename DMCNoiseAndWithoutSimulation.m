load Gz.mat;
D = 2393;
N = 600;
Nu = 1;
lambda =1;
yzad = h2p+0.5;
Tk = 4000;

Fd = [ ones(1,Tk/2).*FDp, ones(1,Tk/2).*(FDp-5) ] ;

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%%Prepare DMC regulator with noise.
dmcnoise = DMCRegWithNoise(Gz,D, N, Nu, lambda);
dmcnoise.reset(F1p, FDp);
dmcnoise.setValue(yzad);

%Prepare DMC regulator.
dmc = DMCReg(Gz,D, N, Nu, lambda);
dmc.reset(F1p);
dmc.setValue(yzad);

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
    uk(k) = dmcnoise.countValue(y(k), Fd(k));
end

figure();
subplot(2,2,1);
plot(ones(Tk,1).*(yzad), 'b-');

hold on;
stairs(y, 'r');
legend('Wyj�cie zadane', 'Wyj�cie regulatora', 'Location', 'southeast');
title(strcat('Regulator DMC z uwzgl�dnieniem zak��ce�'));
xlabel('k');
ylabel('y');
subplot(2,2,3);
stairs( uk, 'g');
hold on;
stairs(Fd, 'r');
legend('Sterowanie', 'Zak��cenie', 'Location', 'east');
xlabel('k');
ylabel('u');
hold off;

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
    uk(k) = dmc.countValue(y(k));
end


subplot(2,2,2);
plot(ones(Tk,1).*(yzad), 'b-');
hold on;
stairs(y, 'r');
legend('Wyj�cie zadane', 'Wyj�cie regulatora', 'Location', 'southeast');
xlabel('k');
ylabel('y');
title(strcat('Regulator DMC konwencjonalny'));
subplot(2,2,4);
stairs( uk, 'g');
hold on;
stairs(Fd, 'r');
legend('Sterowanie', 'Zak��cenie', 'Location', 'east');
xlabel('k');
ylabel('u');
hold off;
print(sprintf('DMCCompareSkok=%.1f.pdf', round(Fd(end)-Fd(1), 1)), '-dpdf');