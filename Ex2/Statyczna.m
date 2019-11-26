y_stat =zeros(2,150);
F1p=73;
FDp = 14;
alfa1 = 20;
alfa2 = 22;

number = 5;

if number == 2
    F1p_array = [(F1p-30) (F1p+30)];
    
    MembershipFunctions = cell(1,size(F1p_array,2));
    MembershipFunctions{1}=MembershipFunction([F1p_array(1)+10 1; F1p_array(2)-10 0]);
    MembershipFunctions{2}=MembershipFunction([F1p_array(1)+10 0; F1p_array(2)-10 1]);
end
if number == 3
    F1p_array = [(F1p-45) F1p (F1p+45)];

    %Prepare membership functions.
    MembershipFunctions = cell(1,size(F1p_array,2));
    MembershipFunctions{1}=MembershipFunction([F1p_array(1)+10 1; F1p_array(2)-10 0]);
    MembershipFunctions{2}=MembershipFunction([F1p_array(1)+10 0; F1p_array(2)-10 1; F1p_array(2)+10 1; F1p_array(3)-10 0]);
    MembershipFunctions{3}=MembershipFunction([F1p_array(2)+10 0; F1p_array(3)-10 1]);
end
if number == 4
    F1p_array = [(F1p-60) (F1p-20) (F1p+20) (F1p+60)];

    %Prepare membership functions.
    MembershipFunctions = cell(1,size(F1p_array,2));
    MembershipFunctions{1}=MembershipFunction([F1p_array(1)+10 1; F1p_array(2)-10 0]);
    MembershipFunctions{2}=MembershipFunction([F1p_array(1)+10 0; F1p_array(2)-10 1; F1p_array(2)+10 1; F1p_array(3)-10 0]);
    MembershipFunctions{3}=MembershipFunction([F1p_array(2)+10 0; F1p_array(3)-10 1; F1p_array(3)+10 1; F1p_array(4)-10 0]);
    MembershipFunctions{4}=MembershipFunction([F1p_array(3)+10 0; F1p_array(4)-10 1]);
end
if number == 5
    F1p_array = [(F1p-60) (F1p-30) F1p (F1p+30) (F1p+60)];

    %Prepare membership functions.
    MembershipFunctions = cell(1,size(F1p_array,2));
    MembershipFunctions{1}=MembershipFunction([F1p_array(1)+10 1; F1p_array(2)-10 0]);
    MembershipFunctions{2}=MembershipFunction([F1p_array(1)+10 0; F1p_array(2)-10 1; F1p_array(2)+10 1; F1p_array(3)-10 0]);
    MembershipFunctions{3}=MembershipFunction([F1p_array(2)+10 0; F1p_array(3)-10 1; F1p_array(3)+10 1; F1p_array(4)-10 0]);
    MembershipFunctions{4}=MembershipFunction([F1p_array(3)+10 0; F1p_array(4)-10 1; F1p_array(4)+10 1; F1p_array(5)-10 0]);
    MembershipFunctions{5}=MembershipFunction([F1p_array(4)+10 0; F1p_array(5)-10 1]);
end
%Main simulation.
%%Prepare DMC regulator.
fuzzy_obj = FuzzyObjectK(F1p_array, MembershipFunctions);

h2_stat = @(x)(x+FDp).^2/(alfa2^2);

for u = 1:150
    [h2, x] = fuzzy_obj.countValue(u); 
    y_stat(2,u) = h2;
end

figure;
plot(1:150,h2_stat(1:150), 1:150, y_stat(2,:));
xlabel('u');
ylabel('y');
title('Porównanie charakterystyki statycznej modelu nieliniowego i rozmytego');
hold on;
scatter(F1p_array, h2_stat(F1p_array), '*r');
legend('Model nieliniowy', 'Model rozmyty', 'Punkty linearyzacji', 'Location', 'southeast');
hold off;
%print(sprintf('pdfs\\ModeleRozmyte\\Porownanie%d.pdf', number), '-dpdf');