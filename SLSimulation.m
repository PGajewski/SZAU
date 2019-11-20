Tk = 4000;
%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
F1p = 73;
h1p=18.9225;
h2p=15.6384;

%Prepare membership functions.
MembershipFunctions = cell(1,5);
MembershipFunctions{1}=MembershipFunction([(h2p_array(1)+2) 1; (h2p_array(2)-2) 0]);
MembershipFunctions{2}=MembershipFunction([(h2p_array(1)+2) 0; (h2p_array(2)-2) 1; (h2p_array(2)+2) 1; (h2p_array(3)-2) 0]);
MembershipFunctions{3}=MembershipFunction([(h2p_array(2)+2) 0; (h2p_array(3)-2) 1; (h2p_array(3)+2) 1; (h2p_array(4)-2) 0]);
MembershipFunctions{4}=MembershipFunction([(h2p_array(3)+2) 0; (h2p_array(4)-2) 1; (h2p_array(4)+2) 1; (h2p_array(5)-2) 0]);
MembershipFunctions{5}=MembershipFunction([(h2p_array(4)+2) 0; (h2p_array(5)-2) 1]);


%%Prepare DMC regulator.
sl = SLReg(2393, 600, 1, 10, -3, 3, ones(600,1)*15, ones(600,1)*20,tfs, MembershipFunctions);
sl.reset(F1p);
sl.setValue(h2p+5);

uk= ones(300,1).*F1p;
y = ones(Tk, 1).*h2p;
h = [h1p, h2p];
%Main simulation loop.
for k=2:Tk
    k
    if k > (300)
        stateHandler = @(t,x) stateFunction(t,x,uk(k - 300), 14);
        [t, h] = ode45(stateHandler,[0 0.5],h(end, :), options);
        y(k) = h(end,2);
    end
    uk(k) = sl.countValue(y(k));
%     if ( k == 2)
%         zmiana = uk(k);
%     end
%     if( uk(k) - zmiana > 10e-6)
%         zmiana = uk(k);
%     end
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