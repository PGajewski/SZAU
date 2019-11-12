function NonlinearPlotsEx1(F1JumpVal, F1JumpTime, FdJumpVal, FdJumpTime, sim_time)
F1p=73;
FDp=14;
tau=150;
h2p=15.6384;
h1p=18.9225;



[t,h] = objectSimulation(F1,FD,tau, sim_time, [h1p, h2p]);

F1vector = ones(sim_time, 1);
F1vector(1:F1(2,1)) = F1(1,2);
F1vector(F1(2,1):sim_time) = F1(2,2);

Fdvector = ones(sim_time, 1);
Fdvector(1:FD(2,1)) = FD(1,2);
Fdvector(FD(2,1):sim_time) = FD(2,2);

subplot(2,1,1);
plot(t,h);
legend('h1','h2 - wyjœcie', 'Location', 'northwest');
xlabel('Czas [s]');
ylabel('Stany wewnetrzne');
title('Przebieg zmiennych stanu');

subplot(2,1,2);
plot(F1vector);
hold on;
plot(Fdvector);
hold off;
legend('F1','Fd - zak³ócenie', 'Location', 'west');
xlabel('Czas [s]');
ylabel('Wejœcia obiektu');
title(sprintf('Przebieg wejœæ'));
print(sprintf('pdfs/SymulacjaObiektu/SymulacjaF1=%dFd=%d.pdf', F1JumpVal, FdJumpVal), '-dpdf');
