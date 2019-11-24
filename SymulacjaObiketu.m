Tk = 1500;


F1p=73;
FDp=14;
tau=150;
h2p=15.6384;
h1p=18.9225;
T = 0.5;

%Prepare ODE options.
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

F1p_array = [(F1p-60) (F1p-30) F1p (F1p+30) (F1p+60)];

%Prepare membership functions.
MembershipFunctions = cell(1,5);
MembershipFunctions{1}=MembershipFunction([F1p_array(1)+10 1; F1p_array(2)-10 0]);
MembershipFunctions{2}=MembershipFunction([F1p_array(1)+10 0; F1p_array(2)-10 1; F1p_array(2)+10 1; F1p_array(3)-10 0]);
MembershipFunctions{3}=MembershipFunction([F1p_array(2)+10 0; F1p_array(3)-10 1; F1p_array(3)+10 1; F1p_array(4)-10 0]);
MembershipFunctions{4}=MembershipFunction([F1p_array(3)+10 0; F1p_array(4)-10 1; F1p_array(4)+10 1; F1p_array(5)-10 0]);
MembershipFunctions{5}=MembershipFunction([F1p_array(4)+10 0; F1p_array(5)-10 1]);


%Main simulation.
%%Prepare DMC regulator.
fuzzy_obj = FuzzyObjectK(F1p_array, MembershipFunctions);

uk= ones(Tk,1).*(F1p+5);
y2 = zeros(Tk, 1);

h1 = [h1p; h2p];
h2 = [h1p; h2p];
%Main simulation loop.

F1=[0, F1p;
    1 F1p+5];
FD=[0, FDp];

[t1,x1,y1] = objectSimulation(F1,FD,tau, Tk*T, h1');

for k=2:Tk    
    %Fuzzy object.
    [temp1, temp2] = fuzzy_obj.countValue(uk(k));
    y2(k) = temp1;
end
figure(8);
plot(t1,y1);
hold on;
plot(T:T:Tk*T,y2);
title('Porownanie obiektu nieliniowego oraz rozmytego');