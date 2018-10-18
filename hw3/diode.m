
clc;
clear all;
close all;

N = 1000;%number of grid points on each region
Nx = 2*N; %total number of grid points
del = (10000)/((N-1)*2+1); %grid spacing such that the discontinuity point
                           %(point where p and n regions touch) not included
delta = del*1e-9;%grid spacing with proper units

%doping profile 
           %units
NA = 2e23; %m-3
ND = 1e23; %m-3
NC = 2.7e25;%m-3
NV=1.04e25;%m-3
%fundamental constants
e0 = 8.85e-12;
epsi = 11.2*e0;
q = 1.6e-19;
vt = 26e-3;
ni = 1e16;

%initialize potential and charge densities
v = zeros(1,Nx);
n = zeros(1,Nx);
p = zeros(1,Nx);
C(1:N) = -1*NA; %C = ND - NA
n(1) = (ni^2)/NA;
p(1) = NA;
p(2:N) = NA/2;
n(2:N) = (ni^2)/(NA/2);
n(Nx) = ND;
p(Nx) = (ni^2)/ND;
n(N+2:Nx) = ND/2;
p(N+2:Nx) = (ni^2)/(ND/2);
C(N+1:Nx) = ND;
rho = q*(p - n + C); %charge densities

emax = 0;
it = 1;

%
% Solve for Equilibruim condition
%
while true
    
    %
    % poisson equation
    %
    for i=2:Nx
        if ((i < Nx))
             vnew(i) = ((v(i-1)+v(i+1))/2) + ((rho(i)*delta^2)/(2*epsi));     
        else
            vnew(i) = v(i-1);
        end
        e = abs((vnew(i)-v(i))/vnew(i));
        if e > emax
            emax = e;
        end
        v(i) = vnew(i);
    end
    
    %
    %Continuity equation using Scharfetter-Gummel scheme
    %
    for i=2:Nx
        if(i<Nx)
            n(i) = n(i-1)*exp((v(i+1)-v(i-1))/(2*vt));
            p(i) = p(i-1)*exp((v(i-1)-v(i+1))/(2*vt));
            rho(i) = q*(p(i) - n(i) + C(i));
        end
    end
    str = sprintf('iteration = %d, max error = %f',it,emax);
    disp(str)
    it=it+1;
    if emax <= 0.01; break; end %tolerance
    emax = 0;
end

Ec = 26e-3*log(NC./n); %eV
Ev = 26e-3*log(p./NV);
Ef(1:Nx) = 0;
%plots
x= 0:delta*1e6:10;

%potential profile
figure
plot(x,v);
xlabel('x (\mum)')
ylabel('Potential (V)')
title('potential profile vs x')
% axis([0 max(x) -1 0.1])

%electron density
figure
plot(x,n.*1e-6);
xlabel('x (\mum)')
ylabel('Electron density (cm^{-3})')
title('Electron density vs x')
axis([0 max(x) -0.5e17 2*max(n)*1e-6])

%hole density
figure
plot(x,p.*1e-6)
xlabel('x (\mum)')
ylabel('Hole density (cm^{-3})')
title('Hole density vs x')
axis([0 max(x) -0.5e17 2*max(p)*1e-6])

%charge density
figure
plot(x(900:1100),rho(900:1100).*1e-6)
xlabel('x (\mum)')
ylabel('charge density (Ccm^{-3})')
title('charge density vs x (zoomed to show the depletion region)')

%band diagram
figure
plot(x(900:1100),Ec(900:1100),x(900:1100),Ev(900:1100),x(900:1100),Ef(900:1100))
legend('Ec','Ef','Ev')
xlabel('x (\mum)')
ylabel('Energy (eV)')
title('Band diagram (zoomed to show the depletion region)')