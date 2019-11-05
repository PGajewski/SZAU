function [output] = MembershipFunction(charPoints, value)
%MembershipFunction Count value of membership function based on their
s = size(charPoints);
for i=2:s(1,1)
    if charPoints(i,1)<=charPoints(i-1,1)
       error("Characteristic points must be rising!");
    end
end

%characteristic points.

%Border conditions.
if charPoints(1,1) < value
    output = charPoints(1,2);
elseif charPoints(end,1) > value
    output = charPoints(end,2);
end 

%Find value
for i=2:s(1,1)
    if charPoints(i,1) < value
        output = (charPoints(i,2)-charPoints(i-1,2))/(charPoints(i,1)-charPoints(i-1,1))*value;
    end
end
end

