y_stat =zeros(2,150);
F1p=73;
FDp = 14;
alfa1 = 20;
alfa2 = 22;
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

h2_stat = @(x)(x+FDp)^2/(alfa2^2);
for u = 1:150
    y_stat(1,u) = h2_stat(u);
    temp = fuzzy_obj.getValue(u); 
    y_stat(2,u) = temp(:,2);
end

figure;
plot(1:150,y_stat(1,:), 1:150, y_stat(2,:));
hold on;
