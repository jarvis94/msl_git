%code to simulate the minority carrier diffusion equation
clc;
clear all;
close all;

%time matrix
T = 0.02e-9; %s
Nt = 500;
t = linspace(0,T,Nt);
dt = t(2) - t(1);
L = 100e-6;

%1D space matrix
Nx= 10001;
x = linspace(-L/2,L/2,Nx);
dx = x(2) - x(1);

%contants
D = 25e-4; %m2/s
k = (D*dt)/dx^2;
% dt = 1.2e-17;
% k = 1/2;

nb = 1e4; % background concentration of minority electrons in P type Si
          %bar
n0 = 1e19; %cm-3 initial dose of electrons at the centre of the bar 

n = nb*ones(length(x),length(t));%concentration matrix of minority carriers

%concentration matrix updated to reflect the initial dose of electrons 
% at the centre spread across 10nm on each side. Choosing 10nm since
% simulation takes a long time if a grid size of 5nm is  
n(x==0,1) = n0; 
i = find(x==0);
n(i+1,1) = n0;
n(i-1,1) = n0;
% n(i+2,1) = 0.5*n0;
% n(i-2,1) = 0.5*n0;

%transformation matrix for implicit finite difference scheme
%A*C(i,t+1) = C(i,t);
%using implicit finite difference scheme since its unconditionaly stable
%i.e no restrictions on dt and dx

A = ((1+2*k)*diag(ones(Nx,1)) - k*diag(ones(Nx-1,1),1) -  k*diag(ones(Nx-1,1),-1));
%Neumann boundary condition
A(1,2) = -2*k;
A(Nx,Nx-1) = -2*k;
Ainv = A^-1;



% solve for concentration matrix for next time instant using
% C(i,t+1) = Ainv*C(i,t);

for j = 1:length(t)-1
    
    n(:,j+1) = Ainv*n(:,j);
%     plot(x(4950:5050)/1e-9,n(4950:5050,j)/n0);
end

%plots zoomed to the centre of the Si bar
figure;
set(gcf,'DefaultLineLineWidth',2);
plot(x(4950:5050)/1e-9,n(4950:5050,1)/n0,x(4950:5050)/1e-9,n(4950:5050,2)/n0,...
     x(4950:5050)/1e-9,n(4950:5050,10)/n0,x(4950:5050)/1e-9,n(4950:5050,200)/n0,...
     x(4950:5050)/1e-9,n(4950:5050,500)/n0);
 axis([-500 500 0 1])
 xlabel('x(nm)')
 ylabel('Normalised concentration')
 tl = [t(1) t(2) t(10) t(200) t(500)]/1e-12;
 str = [repmat('t = ',5,1),num2str(tl'),repmat(' ps',5,1)];
 legend(str)
 title('Normalised Concentration profile of electrons')
 
 %comparison with analytical results
 t1 = 0.36e-12; %time stamp
 M = 30e-9; % total integrated electron dose 

 %constants and prefactor
 Dt = D*t1; 
 pre = M/sqrt(4*pi*Dt);
 
 %analytical solution to 1D diffusion with initial condition being a delta
 %function
 c = pre*exp(-(x.^2)/(4*Dt));
 
 figure
 set(gcf,'DefaultLineLineWidth',2);
 plot(x(4950:5050)/1e-9,n(4950:5050,10)/n0,'o-',x(4950:5050)/1e-9,c(4950:5050))
 xlabel('x(nm)')
 ylabel('Normalised concentration')
 title('Normalised Concentration of electrons at t=0.36ps')
 str1 = sprintf('Numerical');
 str2 = sprintf('Analytical');
 legend(str1,str2)
 axis([-500 500 0 0.3])
 
 