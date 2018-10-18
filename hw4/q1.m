clc;
clear all;
close all;
W = 2e-9;
N = 1000;
x = linspace(0,W,N)';

% Coordinate vector
dx = x(2) - x(1);

%fundamental constants
h = 6.626e-34; %m2kg/s
hbar = h/(2*pi);
m0 = 9.1e-31; %kg
meff = 0.067*m0; %effective mass

%hamiltonion
%discrete double derivative matrix or laplacian matrix
%f''(x) = (f(x-h) - 2f(x) + f(x+h))/h^2
lap = (-2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/(dx^2);

%to enforce boundary conditions
%psi(0),psi(W) = 0 the hamiltonion must be adjusted
lap(1,1) = 0 ; lap(1,2) = 0; lap(2,1) = 0;
lap(N-1,N) = 0; lap(N,N-1) = 0; lap(N,N) = 0;

H=(-hbar^2/(2*meff))*lap; %potential is zero in the well 

%solve the eigen value problem H*psi = E*psi; i.e time independent Schrodinger 
%equation

[psi,E] = eig(H);
Eev = diag(E)/1.6e-19;
figure
plot(x/1e-9,psi(:,3),'r',x/1e-9,psi(:,4),'b'); 
str1 = sprintf('E = %feV',Eev(3));
str2 = sprintf('E = %feV',Eev(4));
xlabel('x (nm)')
ylabel('\psi')
legend(str1,str2)