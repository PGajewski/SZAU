F1p = 73;
number = 2;

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

s = size(MembershipFunctions,2);
figure;
for i=1:s
    v = [];
    F1 = 0:0.1:150;
   for j = F1
       v = [v MembershipFunctions{i}.getValue(j)];
   end
   plot(F1,v);
   hold on;
end
title('Funkcje przynale¿noœci do zbiorów rozmytych');
scatter(F1p_array, ones(1, size(F1p_array,2)));
hold off;
%print(sprintf('pdfs\\ModeleRozmyte\\Funkcja%d.pdf', number), '-dpdf');