clc;
clear all;
close all;
N = 1000;
x1 = linspace(-20,40,N);
x2 = linspace(-20,20,N);
dx1 = x1(2) - x1(1);
dx2 = x2(2) - x2(1);
%linear
u1 = 20*heaviside(-1*x1) + min(heaviside(x1).*x1,20);
%parabolic
u2 = x2.^2;
% plot(x1,u1);
% plot(x2,u2);

%constants
hbar = 1; meff = 1;

%hamiltonion

%discrete double derivative matrix or laplacian matrix
%f''(x) = (f(x-h) - 2f(x) + f(x+h))/h^2
% lap = (-2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/(dx^2);

e = ones(N,1);
lap1 = (spdiags([e -2*e e],[-1 0 1],N,N))/dx1^2; % sparse diagonal form
H1 =(-1*hbar^2/(2*meff))*lap1; %Hamiltonian without potential

lap2 = (spdiags([e -2*e e],[-1 0 1],N,N))/dx2^2; % sparse diagonal form
H2 =(-1*hbar^2/(2*meff))*lap2; %Hamiltonian without potential


%sparse matrix of U
U1 = spdiags(u1',0,N,N);
U2 = spdiags(u2',0,N,N);
%final hamiltonian with potential matrix
H1 = H1 + U1;
H2 = H2 + U2;

%solve for first six Energy states
nmodes = 6;
[psi1,E1] = eigs(H1,nmodes,'sa','Display',0);%'MaxIterations',1000,'Tolerance',1e-10);
[psi2,E2] = eigs(H2,nmodes,'sa','Display',0,'MaxIterations',1000,'Tolerance',1e-10);
E1 = diag(E1);
E2 = diag(E2);

%find the difference between energy states dE(i) = E(i+1) - E(i)
e=ones(6,1);
Diff = spdiags([-1*e e],[0 1],5,6);
dE1 = Diff*E1;
dE2 = Diff*E2;

figure
plot(1:5,dE1,1:5,dE2)
legend('triangular','parabolic')
ylabel('\DeltaE(i)')
xlabel('i')
xt = 1:5;
xticks(xt)
