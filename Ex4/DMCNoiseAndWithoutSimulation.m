load Gz.mat;
D = 2393;
N = 600;
Nu = 1;
lambda =1;
yzad = h2p+2;
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

uk1= ones((Gz.InputDelay(1)),1).*F1p;
y1 = ones(Tk, 1).*h2p;
h = [h1p, h2p];

%Main simulation loop.
for k=2:Tk
    if k > (Gz.InputDelay(1))
        stateHandler = @(t,x) stateFunction(t,x,uk1(k - (Gz.InputDelay(1))), Fd(k));
        [t, h] = ode45(stateHandler,[0 Gz.Ts],h(end, :), options);
        y1(k) = h(end,2);
    end
    uk1(k) = dmcnoise.countValue(y1(k), Fd(k));
end

figure();
subplot(2,2,1);
plot(ones(Tk,1).*(yzad), 'b-');

hold on;
stairs(y1, 'r');
legend('Wyjœcie zadane', 'Wyjœcie regulatora', 'Location', 'southeast');
title(strcat('Regulator DMC z uwzglêdnieniem zak³óceñ'));
xlabel('k');
ylabel('y');
subplot(2,2,3);
stairs( uk1, 'g');
hold on;
stairs(Fd, 'r');
legend('Sterowanie', 'Zak³ócenie', 'Location', 'east');
xlabel('k');
ylabel('u');
hold off;

uk2= ones((Gz.InputDelay(1)),1).*F1p;
y2 = ones(Tk, 1).*h2p;
h = [h1p, h2p];

%Main simulation loop.
for k=2:Tk
    if k > (Gz.InputDelay(1))
        stateHandler = @(t,x) stateFunction(t,x,uk2(k - (Gz.InputDelay(1))), Fd(k));
        [t, h] = ode45(stateHandler,[0 Gz.Ts],h(end, :), options);
        y2(k) = h(end,2);
    end
    uk2(k) = dmc.countValue(y2(k));
end


subplot(2,2,2);
plot(ones(Tk,1).*(yzad), 'b-');
hold on;
stairs(y2, 'r');
legend('Wyjœcie zadane', 'Wyjœcie regulatora', 'Location', 'southeast');
xlabel('k');
ylabel('y');
title(strcat('Regulator DMC konwencjonalny'));
subplot(2,2,4);
stairs( uk2, 'g');
hold on;
stairs(Fd, 'r');
legend('Sterowanie', 'Zak³ócenie', 'Location', 'east');
xlabel('k');
ylabel('u');
hold off;
%print(sprintf('DMCCompareSkok=%.1f.pdf', round(Fd(end)-Fd(1), 1)), '-dpdf');

figure;
subplot(2,1,1)
plot(ones(Tk,1).*(yzad), 'b-');
hold on;
stairs(y2, 'r');
stairs(y1, 'g');
legend('Wyjœcie zadane', 'Wyjœcie regulatora bez', 'Wyjœcie regulatora z', 'Location', 'southeast');
title(strcat('Porównanie regulatorów DMC z i bez uwzglêdniania zak³óceñ'));
xlabel('k');
ylabel('y');
subplot(2,1,2);
stairs( uk2, 'r');
hold on;
stairs( uk1, 'g');
stairs(Fd, 'b');
legend('Sterowanie bez', 'Sterowanie z', 'Zak³ócenie', 'Location', 'east');
xlabel('k');
ylabel('u');
hold off;