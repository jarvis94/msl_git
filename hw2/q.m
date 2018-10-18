clc;
clear all;
close all;
l = 700;
h = 1;
hr= h*1e-9;

NAm = 1e23;
NDp = 5e23;
xp = 500;
xn = 200;
nx = (l/h)+1;
nxmid = 251;
nxp =xp/h+1;
nxn = xn/h+1;
emax = 0;
v=zeros(1,nx);
e0 = 8.85e-12;
epsi = 11.2*e0;
q = 1.6e-19;

%tri diagonal matrix
A = zeros(nx,nx);
A(1,1) = 1;
for i=2:nx-1
  A(i,i) = -2;
  if (i>2)
      A(i,i-1) = 1;
  end
  A(i,i+1) = 1;
end
A(nx,nx) = -1;
A(nx,nx-1) = 1;

%charge density matrix

C =zeros(1,nx);
C(1:nxn-1)=0;
C(nxn:nxmid-1) = -1*((q*NDp*hr^2)/(epsi));
C(nxmid:nxp-1)= ((q*NAm*hr^2)/(epsi));
C(nxp:nx) = 0;
V=(A^-1)*C';

subplot(3,1,1)

rho = C.*((-1*epsi)/(q*hr^2));
plot(1:nx,rho);
title('uniform doping')
ylabel('charge density (C/m^{-3})')
axis([0 800 -2e23 6e23]);
subplot(3,1,2)
plot(1:nx,V);
ylabel('potential (V)')
E = -1*gradient(V);
subplot(3,1,3)
plot(1:nx,E);
ylabel('Electric Field (V/nm)')
xlabel('position (nm)')


%linear grading
mn = 1e22;
mp = 4e20;
n = 0:50;
fn = -1*((q*mn.*n*hr^2)/(epsi));
p = 250:-1:0;
fp = ((q*mp.*p*hr^2)/(epsi));

C(nxn:nxmid-1) = fn(1:50);
trapz(1:50,fn(1:50))
C(nxmid:nxp-1)= fp(2:251);
trapz(2:251,fp(2:251))
V2=(A^-1)*C';
figure
subplot(3,1,1)
rho = C.*((-1*epsi)/(q*hr^2));
plot(1:nx,rho);
title('Linearly Graded Junction')
axis([0 800 -2e23 6e23]);
ylabel('charge density (C/m^{-3})')
subplot(3,1,2)
plot(1:nx,V2);
ylabel('potential (V)')
subplot(3,1,3)
E2 = -1*gradient(V2);
plot(1:nx,E2);
ylabel('Electric Field (V/nm)')
xlabel('position (nm)')

%gaussian

x=0:50;
x1 = 0:250;
norm = normpdf(x,25,7);
trapz(x(1:50),norm(1:50))
norm1 = normpdf(x1,125,40);
trapz(x1(2:250),norm1(2:250))

fn = (-1*((q*NDp*50*hr^2)/(epsi))).*norm;
fp = (((q*NAm*250*hr^2)/(epsi))).*norm1;

C(nxn) = 0;
C(nxn+1:nxmid-1) = fn(2:50);
C(nxmid) = 0;
C(nxmid+1:nxp-1)= fp(2:250);

V3=(A^-1)*C';
figure
subplot(3,1,1)
rho = C.*((-1*epsi)/(q*hr^2));
plot(1:nx,rho);
title('Gaussian Doping profile')
ylabel('charge density (C/m^{-3})')
subplot(3,1,2)
plot(1:nx,V3);
ylabel('potential (V)')
E3 = -1*gradient(V3);
subplot(3,1,3)
plot(1:nx,E3);
ylabel('Electric Field (V/nm)')
xlabel('position (nm)')

