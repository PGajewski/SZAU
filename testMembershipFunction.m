MembershipFunctions = cell(1,5);
MembershipFunctions{1}=MembershipFunction([(F1p-50) 1; (F1p-40) 0]);
MembershipFunctions{2}=MembershipFunction([(F1p-50) 0; (F1p-40) 1; (F1p-20) 1; (F1p-10) 0]);
MembershipFunctions{3}=MembershipFunction([(F1p-20) 0; (F1p-10) 1; (F1p+10) 1; (F1p+20) 0]);
MembershipFunctions{4}=MembershipFunction([(F1p+10) 0; (F1p+20) 1; (F1p+40) 1; (F1p+50) 0]);
MembershipFunctions{5}=MembershipFunction([(F1p+40) 0; (F1p+50) 1]);

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