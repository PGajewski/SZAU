%% Wyznaczanie charakterystyki statycznej obiektu nieliniowego.
y_stat =zeros(2,200);
sim_time = 2000;
Ts = 0.5;

h2_stat = @(x)(x+FDp)^2/(alfa2^2)
for u = 1:150
    u
    y_stat(1,u) = h2_stat(u);
    
    fuzzy_obj.reset([0.1;0.1],[0.1;0.1]);
    temp = 0;
    for t=Ts:Ts:sim_time/Ts
       temp = fuzzy_obj.countValue([u;FDp]); 
    end
    y_stat(2,u) = temp;
end

figure;
plot(1:200,y_stat(1,:), 1:200, y_stat(2,:));
hold on;
