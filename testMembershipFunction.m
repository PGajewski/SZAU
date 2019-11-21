MembershipFunctions = cell(1,5);
MembershipFunctions{1}=MembershipFunction([F1p_array(1)+10 1; F1p_array(2)-10 0]);
MembershipFunctions{2}=MembershipFunction([F1p_array(1)+10 0; F1p_array(2)-10 1; F1p_array(2)+10 1; F1p_array(3)-10 0]);
MembershipFunctions{3}=MembershipFunction([F1p_array(2)+10 0; F1p_array(3)-10 1; F1p_array(3)+10 1; F1p_array(4)-10 0]);
MembershipFunctions{4}=MembershipFunction([F1p_array(3)+10 0; F1p_array(4)-10 1; F1p_array(4)+10 1; F1p_array(5)-10 0]);
MembershipFunctions{5}=MembershipFunction([F1p_array(4)+10 0; F1p_array(5)-10 1]);
s = size(MembershipFunctions);
figure;
for i=1:5
    v = [];
    F1 = 0:0.1:150;
   for j = F1
       v = [v MembershipFunctions{i}.getValue(j)];
   end
   plot(F1,v);
   hold on;
end
hold off;