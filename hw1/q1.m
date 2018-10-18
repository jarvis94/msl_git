clc;
clear all;
l = 1000;
d = 10;
[cap,cap_ideal,par_cap]=cap_solver(l,d);
s1 = sprintf("Net capacitance: %d",cap);
s2 = sprintf("Parasictic cap: %d",par_cap);
disp(s1)
disp(s2)