function [y, u] = DMC(Gz, D, N, Nu, lambda, hlin, ulin, Tk, yzad)
% Gz - transmitancja dyskretna
% D - hotyzont dynamiki, (1200)
% N - horyzont predykcji,
% Nu - d³ugoœæ horyzontu sterownia
% lambda - parametr regulatora (1200)
% hlin - 2x1 lub 1x2 wektor wartoœci h pracy
% ulin - -||- u pracy
% Tk - d³ugoœæ symulacji
% yzad - wartoœæ zadana


u(1:(Gz.InputDelay(1))) = 0;
y1(1:(Gz.InputDelay(1))) = 0;
u((Gz.InputDelay(1)):((Gz.InputDelay(1))+D)) = 1;
s = zeros(D, 1);


% Wyznaczanie wektora s
for k=(Gz.InputDelay(1))+1:(Gz.InputDelay(1))+D
    y1(k)=-Gz.Denominator{1}(2)*y1(k-1)-Gz.Denominator{1}(3)*y1(k-2)+Gz.Numerator{1}(2)*u(k-((Gz.InputDelay(1))-1))+Gz.Numerator{1}(3)*u(k-(Gz.InputDelay(1)));
    s(k-(Gz.InputDelay(1)))=y1(k);
end

psi = eye(N);
LAMBDA = eye(Nu)*lambda;


% Wyznaczanie macierzy M
M = zeros(N, Nu);
for i=1:1:Nu
    M(i:N,i)=s(1:N-i+1)';
end


% Wyznaczanie macierzy Mp
Mp= zeros(N, D-1);
for i=1:N
   for j=1:D-1
      if i+j<=D
         Mp(i,j)=s(i+j)-s(j);
      else
         Mp(i,j)=s(D)-s(j);
      end   
   end
end

% Wektor wzmocnieñ - wyznaczany raz (offline)
K=(M'*psi*M+LAMBDA)^(-1)*M'*psi;

% Przygotowanie zmiennych do symulacji
yzad(1:N)= yzad;
yzad=yzad';
deltaup=zeros(1,D-1)';
uk= ones((Gz.InputDelay(1)),1).*ulin(1);
y = ones(Tk, 1).*hlin(2);
h1 = hlin(1);


% g³ówna pêtla symulacji
for k=2:Tk
    if k > (Gz.InputDelay(1))
        h = stateFunction(0,[h1, y(k-1)],uk(k - (Gz.InputDelay(1))), ulin(2));
        h1 = h1 + h(1);
        y(k) = y(k-1) + h(2);
    end
    
    
    % aktualizacja wektora aktualnej wartoœci wyjœcia
    yk=ones(N,1)*y(k);
    
    % wyliczenie nowego wektora odpowiedzi swobodnej
    y0=yk+Mp*deltaup;
    
    % wyliczenie wektora zmian sterowania
    deltauk=K*(yzad-y0);
    
    % prawo regulacji (w chwili 1. uk(0) = 0)
    if k==1
        uk(k)=deltauk(1);
    else
        uk(k)=uk(k-1)+deltauk(1);
    end
    
    % aktualizacja poprzednich zmian sterowania
    deltaup=[deltauk(1) deltaup(1:end-1)']';
end

figure;
plot( uk, 'g');
hold on;
stairs(yzad, 'b');
stairs(y, 'r');
title(strcat('Dzialanie regulatora dla nastaw D=', num2str(D), ' N =', num2str(N), ' Nu=', num2str(Nu), ' lambda=',num2str(lambda)));
grid on;
legend('Sterowanie', 'Wyjœcie zadane', 'wyjœcie regulatora', 'Location', 'east');
xlabel('k');
ylabel('y/u');