clc;
clear all;
close all;
l1= 20e-9;
l2= 5e-9;
l3= 20e-9;
L = l1+l2+l3;
N=1000;
x = linspace(0,L,N);
dx = x(2) - x(1);
%constants
h = 6.62607004e-34; %m2kg/s
hbar = h/(2*pi);
m0 = 9.1e-31; %kg
me_GaAs = 0.067*m0; %effective mass electron in GAs
me_AlAs = 0.150*m0; %effective mass electrom in AlAs
k = 0.6;
me_inv = (k/me_AlAs) + ((1-k)/me_GaAs);
meff = 1/me_inv; %effective mass of electron AlxGa1-xAs

%Potentials
Eg_1 = 1.41; % bandgap of GaAs eV
Eg_2 = 1.908; % bandgap of AlxGa1-xAs
dEg = Eg_2 - Eg_1;
dEc = 0.6*dEg;
dEv = 0.4*dEg;

%potentials are defined considering the conduction band of GaAs as
%reference

Uc = dEc*(heaviside(l1-x)+heaviside(x-l1-l2));
Uv = -dEv*(heaviside(l1-x)+heaviside(x-l1-l2))-Eg_1*ones(1,N);
figure
plot(x,Uc,x,Uv);


%hamiltonion
%discrete double derivative matrix or laplacian matrix
%f''(x) = (f(x-h) - 2f(x) + f(x+h))/h^2
% lap = (-2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/(dx^2);
e = ones(N,1);
lap = (spdiags([e -2*e e],[-1 0 1],N,N))/dx^2; % sparse diagonal form

% why are not enforcing boundary conditions here?
%to enforce boundary conditions
%psi(0),psi(W) = 0 the hamiltonion must be adjusted
% lap(1,1) = 0 ; lap(1,2) = 0; lap(2,1) = 0;
% lap(N-1,N) = 0; lap(N,N-1) = 0; lap(N,N) = 0;

H=(-1*hbar^2/(2*meff))*lap; %Hamiltonian without potential

% Ucs = 1.6e-19*diag(Uc);
% Ucs(1,1) = 0; Ucs(N,N) = 0;
%sparse matrix of Uc
Uvs = 1.6e-19*spdiags(-1*Uv',0,N,N);%in Joules
%final hamiltonian with potential
H = H + Uvs;

nmodes = 1;
[psi,E] = eigs(H,nmodes,'sa','Display',0);

Usc = 1.6e-19*Uv*max(abs(psi(:)))/max(abs(1.6e-19*Uv));
figure
plot(x,psi,x,Usc,'--');
E_ev = E(1)/1.6e-19;
str = sprintf('E = %f',E_ev);
legend(str)
figure
plot(Usc)