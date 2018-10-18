clc;
clear all;
l = 10:50:1000;
stop = size(l);
d = 10;
for i=1:stop(2)
    disp(l(i))
    [cap(i),cap_ideal(i),par_cap(i)]=cap_solver(l(i),d);
    
end
plot(l,cap)
xlabel('length of the plates')
ylabel('Net capacitance')
title('capacitance vs length')
figure
plot(l,par_cap);
xlabel('length of the plates')
ylabel('Net capacitance')
title('parasitic capacitance vs length')
figure
plot(l,(par_cap./cap).*100)
xlabel('length of the plates')
ylabel('% parasitic capacitance')
title('% parasitic capacitance vs length')