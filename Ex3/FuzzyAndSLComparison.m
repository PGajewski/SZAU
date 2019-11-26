%% Fuzzy DMC.
load Gz.mat;
Tk = 5000;
yzad = h2p+8;
Dd = 2393;
N = 600;
Nu = 1;
lambda =1;
%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

tau= 150;
C1 = 0.35;
C2 = 0.3;
alfa1 = 20;
alfa2 = 22;

F1p_array = [(F1p-60) (F1p-30) F1p (F1p+30) (F1p+60)];
h1p_array = [];
h2p_array = [];
tfs = cell(size(F1p_array));
sss = cell(size(F1p_array));
T=0.5;



i = 1;
for local_F1p=F1p_array
    %Simulate non linear object to get h vector of new work point.
    x_p = getLinearModel(local_F1p,FDp);
    h1p_array = [h1p_array x_p(1)];
    h2p_array = [h2p_array x_p(2)];
    %Linearization.
    a11 = alfa1 * x_p(1)^(-5/2)/(2*C1) - 2/(3*C1) * (local_F1p + FDp) * x_p(1)^(-3);
    a12 = 0;
    a21 = alfa1/(6*C2) * x_p(1)^(-1/2) * x_p(2)^-2;
    a22 = -2 * alfa1/(3*C2) * sqrt(x_p(1)) * x_p(2)^-3 + alfa2/(2*C2) * x_p(2)^(-5/2);
    u11 = x_p(1)^(-2)/(3*C1);
    u12 = x_p(1)^(-2)/(3*C1);
    u21 = 0;
    u22 = 0;
    
    A = [a11, a12; a21, a22];
    B = [u11 u12; u21 u22];
    C = [0 1];
    D = [0 0];

    tfs{i} = c2d(tf(ss(A,B,C,D,'InputDelay',[tau,0])),T,'zoh');
    sss{i} = c2d(ss(A,B,C,D),T,'zoh');
    sss{i}.InputDelay = [tau/T,0];
    i = i + 1;
end
D = Dd;
%Prepare local DMCs.
LocalDMCs = cell(1,5);
LocalDMCs{1}=DMCReg(tfs{1},2393, 600, 1, 1);
LocalDMCs{2}=DMCReg(tfs{2},D, N, Nu, lambda);
LocalDMCs{3}=DMCReg(tfs{3},D, N, Nu, lambda);
LocalDMCs{4}=DMCReg(tfs{4},D, N, Nu, lambda);
LocalDMCs{5}=DMCReg(tfs{5},D, N, Nu, lambda);

%Prepare membership functions.
MembershipFunctions = cell(1,5);
MembershipFunctions{1}=MembershipFunction([(h2p_array(1)+2) 1; (h2p_array(2)-2) 0]);
MembershipFunctions{2}=MembershipFunction([(h2p_array(1)+2) 0; (h2p_array(2)-2) 1; (h2p_array(2)+2) 1; (h2p_array(3)-2) 0]);
MembershipFunctions{3}=MembershipFunction([(h2p_array(2)+2) 0; (h2p_array(3)-2) 1; (h2p_array(3)+2) 1; (h2p_array(4)-2) 0]);
MembershipFunctions{4}=MembershipFunction([(h2p_array(3)+2) 0; (h2p_array(4)-2) 1; (h2p_array(4)+2) 1; (h2p_array(5)-2) 0]);
MembershipFunctions{5}=MembershipFunction([(h2p_array(4)+2) 0; (h2p_array(5)-2) 1]);

%Main simulation.
%%Prepare DMC regulator.
fuzzy_dmc = FuzzyDMCReg(LocalDMCs, MembershipFunctions);
fuzzy_dmc.reset(F1p);
fuzzy_dmc.setValue(yzad);

uk1= ones((LocalDMCs{1}.Gz.InputDelay(1)),1).*F1p;
y1 = ones(Tk, 1).*h2p;
h = [h1p, h2p];


%Main simulation loop.
for k=2:Tk
    if k > (LocalDMCs{1}.Gz.InputDelay(1))
        stateHandler = @(t,x) stateFunction(t,x,uk1(k - (LocalDMCs{1}.Gz.InputDelay(1))), FDp);
        [t, h] = ode45(stateHandler,[0 LocalDMCs{1}.Gz.Ts],h(end, :), options);
        y1(k) = h(end,2);
    end
    uk1(k) = fuzzy_dmc.countValue(y1(k));
end


%%Prepare DMC regulator.
sl = SLReg(2393, 600, 1, 10, -4, 4, ones(600,1)*10, ones(600,1)*25,tfs, MembershipFunctions);
sl.reset(F1p);
sl.setValue(yzad);

uk2= ones(300,1).*F1p;
y2 = ones(Tk, 1).*h2p;
h = [h1p, h2p];
%Main simulation loop.
for k=2:Tk
    k
    if k > (300)
        stateHandler = @(t,x) stateFunction(t,x,uk2(k - 300), 14);
        [t, h] = ode45(stateHandler,[0 0.5],h(end, :), options);
        y2(k) = h(end,2);
    end
    uk2(k) = sl.countValue(y2(k));
end


figure();
subplot(2,1,1);
stairs(ones(Tk,1).*(yzad), 'g');
hold on;
stairs(y1, 'r');
stairs(y2, 'b');
title('Porównanie rozmytego i SL dla nastaw D=2393, N =600, Nu=1, lambda=1');
legend('Wyjœcie zadane', 'Wyjœcie regulatora roz.', 'Wyjœcie regulatora SL' , 'Location', 'east');
xlabel('k');
ylabel('y');
subplot(2,1,2);
stairs( uk1, 'r');
hold on;
stairs( uk2, 'b');
xlabel('k');
ylabel('u');
legend('Sterowanie roz.', 'Sterowanie SL');
hold off;