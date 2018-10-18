clc;
clear all;
close all;

T = 200e-9; %s
Nt = 100;
t = linspace(0,T,Nt);
dt = t(2) - t(1);

L = 100e-6; 
Nx= 1001;
x = linspace(-L/2,L/2,Nx);
dx = x(2) - x(1);
 
D = 25e-4; %m2/s
k = (D*dt)/dx^2;

n = zeros(length(x),length(t));
n(x==0,1) = 1e13; %m-3/s 
A = ((1-2*k)*diag(ones(Nx,1)) + k*diag(ones(Nx-1,1),1) +  k*diag(ones(Nx-1,1),-1));
A(1,1) = 0 ; A(1,2) = 0; A(2,1) = 0;
A(Nx-1,Nx) = 0; A(Nx,Nx-1) = 0; A(Nx,Nx) = 0;


for j = 2:length(t)
    n(:,j+1) = A*n(:,j);
end