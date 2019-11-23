T=0.5;
Tk = 1500;

F1p_array = [(F1p-60) (F1p-30) F1p (F1p+30) (F1p+60)];
h1p_array = [];
h2p_array = [];
sim_time = 2000;
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
h2_stat = @(x)(x+FDp).^2/(alfa2^2);

% plot(F1p-60:0.1:F1p+60, h2_stat(F1p-60:0.1:F1p+60));
% hold on;
% scatter(F1p_array, h2p_array);
%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%Prepare local DMCs.
LocalObjects = cell(1,5);
LocalObjects{1}=sss{1};
LocalObjects{2}=sss{2};
LocalObjects{3}=sss{3};
LocalObjects{4}=sss{4};
LocalObjects{5}=sss{5};

%Prepare membership functions.
MembershipFunctions = cell(1,5);
MembershipFunctions{1}=MembershipFunction([F1p_array(1)+10 1; F1p_array(2)-10 0]);
MembershipFunctions{2}=MembershipFunction([F1p_array(1)+10 0; F1p_array(2)-10 1; F1p_array(2)+10 1; F1p_array(3)-10 0]);
MembershipFunctions{3}=MembershipFunction([F1p_array(2)+10 0; F1p_array(3)-10 1; F1p_array(3)+10 1; F1p_array(4)-10 0]);
MembershipFunctions{4}=MembershipFunction([F1p_array(3)+10 0; F1p_array(4)-10 1; F1p_array(4)+10 1; F1p_array(5)-10 0]);
MembershipFunctions{5}=MembershipFunction([F1p_array(4)+10 0; F1p_array(5)-10 1]);

%Prepare x_0 of models.
x_0_models = cell(1,5);
x_0_models{1} = [h1p_array(1); h2p_array(1)];
x_0_models{2} = [h1p_array(2); h2p_array(2)];
x_0_models{3} = [h1p_array(3); h2p_array(3)];
x_0_models{4} = [h1p_array(4); h2p_array(4)];
x_0_models{5} = [h1p_array(5); h2p_array(5)];

%Prepare x_0 of models.
u_0_models = cell(1,5);
u_0_models{1} = [F1p_array(1); FDp];
u_0_models{2} = [F1p_array(2); FDp];
u_0_models{3} = [F1p_array(3); FDp];
u_0_models{4} = [F1p_array(4); FDp]; 
u_0_models{5} = [F1p_array(5); FDp];

%Main simulation.
%%Prepare DMC regulator.
fuzzy_obj = FuzzyObject(LocalObjects, MembershipFunctions,x_0_models, u_0_models);
fuzzy_obj.reset([h1p;h2p],[F1p;FDp]);

uk= ones(Tk,1).*(F1p+50);
y2 = zeros(Tk, 1);

h1 = [h1p; h2p];
h2 = [h1p; h2p];
%Main simulation loop.

F1=[0, F1p;
    1 F1p+50];
FD=[0, FDp];

[t1,x1,y1] = objectSimulation(F1,FD,tau, Tk*T, h1');

for k=2:Tk    
    %Fuzzy object.
    [temp1, temp2] = fuzzy_obj.countValue([uk(k);FDp]);
    y2(k) = temp1;
    h2 = [h2 temp2];
    k
end
figure(8);
plot(t1,y1);
hold on;
plot(T:T:Tk*T,y2);
title('Porownanie obiektu nieliniowego oraz rozmytego');