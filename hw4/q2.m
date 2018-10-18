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
h = 6.626e-34; %m2kg/s
hbar = h/(2*pi);
m0 = 9.1e-31; %kg
me_GaAs = 0.067*m0; %effective mass electron in GAs
me_AlAs = 0.150*m0; %effective mass electrom in AlAs

mh_GaAs = 0.45*m0; %effective mass electron in GAs
mh_AlAs = 0.76*m0; %effective mass electrom in AlAs

k = 0.6;
me_inv = (k/me_AlAs) + ((1-k)/me_GaAs);
meff = 1/me_inv; %effective mass of electron AlxGa1-xAs

mh_inv = (k/me_AlAs) + ((1-k)/me_GaAs);
mhff = 1/me_inv; %effective mass of electron AlxGa1-xAs

%Potentials
Eg_1 = 1.41; % bandgap of GaAs eV
Eg_2 = 1.908; % bandgap of AlxGa1-xAs
dEg = Eg_2 - Eg_1;
dEc = 0.6*dEg;
dEv = 0.4*dEg;

%potentials are defined considering the conduction band of GaAs as
%reference

Uc = dEc*(heaviside(l1-x)+heaviside(x-l1-l2));
Uv = dEv*(heaviside(l1-x)+heaviside(x-l1-l2));
% figure
% plot(x,Uc,x,Uv);


%hamiltonion
%discrete double derivative matrix or laplacian matrix
%f''(x) = (f(x-h) - 2f(x) + f(x+h))/h^2
%  lap = (-2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/(dx^2);
e = ones(N,1);
lap = (spdiags([e -2*e e],[-1 0 1],N,N))/dx^2; % sparse diagonal form


%to enforce boundary conditions
%psi(0),psi(W) = 0 the hamiltonion must be adjusted
% lap(1,1) = 0 ; lap(1,2) = 0; lap(2,1) = 0;
% lap(N-1,N) = 0; lap(N,N-1) = 0; lap(N,N) = 0;


%prefactor for electrons
hn1 = ones(N,1);
hn1(1:446) = (-1*hbar^2/(2*meff));
hn1(447:555) = (-1*hbar^2/(2*me_GaAs));
hn1(556:end) = (-1*hbar^2/(2*meff));

%prefactor for holes
hp1 = ones(N,1);
hp1(1:446) = (-1*hbar^2/(2*mhff));
hp1(447:555) = (-1*hbar^2/(2*mh_GaAs));
hp1(556:end) = (-1*hbar^2/(2*mhff));

Hn=hn1.*lap; %Hamiltonian without potential electron

Hp=hp1.*lap; %Hamiltonian without potential holes


%electron potential well
Ucs = 1.6e-19*diag(Uc);

%hole potential well
Uvs = 1.6e-19*diag(Uv);

% Ucs(1,1) = 0; Ucs(N,N) = 0;
%sparse matrix of Uc
%Ucs = 1.6e-19*spdiags(Uc',0,N,N);%in Joules
%final hamiltonian with potential
Hn = Hn + Ucs;
Hp = Hp + Uvs;

nmodes = 3;

%solve for electroms
[psie,Ee] = eigs(Hn,nmodes,'sa','Display',0,'Maxiterations',500,'Tolerance',1e-10); 
[Ee,inde] = sort(diag(Ee));


%solve for holes
[psip,Ep] = eigs(Hp,nmodes,'sa','Display',0,'Maxiterations',500,'Tolerance',1e-10); 
[Ep,indp] = sort(diag(Ep));

Usc = 1.6e-19*Uc*max(abs(psie(:)))/max(abs(1.6e-19*Uc));
Usv = 1.6e-19*Uv*max(abs(psip(:)))/max(abs(1.6e-19*Uv));

figure
plot(x/1e-9,-1*psie(:,1),x/1e-9,Usc,'--');
En_ev = Ee(1)/1.6e-19;
str = sprintf('E = %f',En_ev);
legend(str)
xlabel('x (nm)')
ylabel('\psi')
title('Electron wave function')
figure
plot(x/1e-9,psip(:,1),x/1e-9,Usv,'--');
Ep_ev = Ep(1)/1.6e-19;
str = sprintf('E = %f',Ep_ev);
legend(str)
xlabel('x (nm)')
ylabel('\psi')
title('Hole wave function')
figure

%plotting energy levels
Ee_ev = Ee/1.6e-19;
Ep_ev = Ep/1.6e-19;
Uv = -1*Uv - Eg_1*ones(1,N);

Ep_ev = -1*Ep_ev - Eg_1;
plot(x,Uc,'k',x(446:555),Ee_ev*ones(1,110),'--',x(446:555),Ep_ev*ones(1,110),'--',x,Uv,'k');
lgnd_str1 = [repmat('E = ',nmodes,1),num2str(Ee_ev)];
lgnd_str2 = [repmat('E = ',nmodes,1),num2str(Ep_ev)];
legend(lgnd_str1,lgnd_str2)
xlabel('x')
ylabel('Energy (eV)')
title('Bang diagram with quantized energy levels in the well')

%density of states and carrier profile
Ef = 0;
kt = 26e-3;
E = 0:0.01:2*dEc;
for i = 1:length(E)
    f(i) = 1/(exp((E(i)-Ef)/kt)+1);
end

NE = zeros(1,length(E));
for i = 1:nmodes
    NE(1:length(E)) = NE + (me_GaAs/(hbar^2*pi))*heaviside(E-Ee_ev(i));
end
n = NE.*f;
figure
plot(E,f)
xlabel('Energy(eV)');
ylabel('Probablity');
figure
plot(E,NE)
xlabel('Energy(eV)');
ylabel('Density of states');
figure
plot(E,n)
% line([Ee_ev Ee_ev],[min(n) max(n)])
xlabel('Energy (eV)')
ylabel('Electron Density n(E)')
Ep = -1*Eg_1:-0.01:-2*dEv-Eg_1;
for i = 1:length(Ep)
    fp(i) = (1/(exp((Ef-Ep(i))/kt)+1));
end

NEp = zeros(1,length(Ep));
for i = 1:nmodes
    NEp(1:length(Ep)) = NEp + (mh_GaAs/(hbar^2*pi))*heaviside(Ep_ev(i)-Ep);
end
p = NEp.*fp;
figure
plot(Ep,fp)
xlabel('Energy(eV)');
ylabel('Probablity');
figure
plot(Ep,NEp)
xlabel('Energy(eV)');
ylabel('Density of states');
figure
plot(Ep,p)
% line([Ep_ev Ep_ev],[min(p) max(p)])
xlabel('Energy (eV)')
ylabel('Electron Density n(E)')