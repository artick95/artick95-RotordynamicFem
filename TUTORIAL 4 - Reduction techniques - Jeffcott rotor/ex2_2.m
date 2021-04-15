clear all
clc

E=2e11;
L=.3;
d=2.5e-3;
m=.5;
Jt=.01;
Jp=.02;
g=-9.81;
Fa=g*m;
I=pi*d^4/64;
delta=(Jp-Jt)/m;

% Jp>Jt: disk rotor and only one critical speed exists
% r/l = 0.0042 hence consider Euler-Bernoulli beam (phi=0)
K=E*I/L^3*[12   6*L     -12     6*L
           6*L  4*L^2   -6*L    2*L^2
           -12  -6*L    12      -6*L
           6*L  2*L^2   -6*L    4*L^2];
       
Kg=-Fa/(30*L)*[36    3*L     -36     3*L
              3*L   4*L^2   -3*L    -L^2
              -36   -3*L    36      -3*L
              3*L   -L^2    -3*L    4*L^2];
       
% M=[m    0   0   0
%    0    Jt  0   0
%    0    0   m   0
%    0    0   0   Jt];
% 
% G=[0    0   0   0
%    0    0   0   -Jp
%    0    0   0   0
%    0    Jp  0   0];

% Mq''-i*omega*G*q'+Ktot*q=0

Kxy=K;
Kxy(1:2,3:4)=zeros(2);
Kxy(3:4,1:2)=zeros(2);
Kgxy=Kg;
Kgxy(1:2,3:4)=zeros(2);
Kgxy(3:4,1:2)=zeros(2);


%Plotting the Campbell diagram
Omega=linspace(0,60,61);
wT=zeros(4,size(Omega,2));
wC=wT;
wN=wC;

% Mass below the contraint (tension)
Ktot=Kxy+Kgxy;         
wcrT1=1/sqrt(2*m)*sqrt(Ktot(1,1)-Ktot(2,2)/delta+sqrt((Ktot(1,1)+Ktot(2,2)/delta)^2-4*Ktot(1,2)^2/delta));
wcrT2=1/sqrt(2*m)*sqrt(Ktot(1,1)-Ktot(2,2)/delta-sqrt((Ktot(1,1)+Ktot(2,2)/delta)^2-4*Ktot(1,2)^2/delta)); %complex

disp('Critical speed under tensile stress [rad/s] = ')
disp(wcrT1)
disp(' ')

for j=1:size(Omega,2)
    wT(:,j)=roots([1 -Omega(j)*Jp/Jt -(Ktot(1,1)/m+Ktot(2,2)/Jt) +Omega(j)*Ktot(1,1)*Jp/(m*Jt) (Ktot(1,1)*Ktot(2,2)-Ktot(1,2)^2)/(m*Jt)]);
end
Omega_star=Omega/wcrT1;
wT_star=wT/wcrT1;
plot(Omega_star,real(wT_star),Omega_star,Omega_star,1,Omega_star,Omega_star,Jp/Jt*Omega_star)

% Mass above the contraint (compression)
Kgxy=-Kgxy;

Ktot=Kxy+Kgxy;          
wcrC1=1/sqrt(2*m)*sqrt(Ktot(1,1)-Ktot(2,2)/delta+sqrt((Ktot(1,1)+Ktot(2,2)/delta)^2-4*Ktot(1,2)^2/delta));
wcrC2=1/sqrt(2*m)*sqrt(Ktot(1,1)-Ktot(2,2)/delta-sqrt((Ktot(1,1)+Ktot(2,2)/delta)^2-4*Ktot(1,2)^2/delta)); %complex

disp('Critical speed under compressive stress [rad/s] = ')
disp(wcrC1)
disp(' ')

figure 
for j=1:size(Omega,2)
    wC(:,j)=roots([1 -Omega(j)*Jp/Jt -(Ktot(1,1)/m+Ktot(2,2)/Jt) +Omega(j)*Ktot(1,1)*Jp/(m*Jt) (Ktot(1,1)*Ktot(2,2)-Ktot(1,2)^2)/(m*Jt)]);
end
Omega_star=Omega/wcrC1;
wC_star=wC/wcrC1;
plot(Omega_star,real(wC_star),Omega_star,Omega_star,1,Omega_star,Omega_star,Jp/Jt*Omega_star)

%No gravity
Ktot=Kxy;
wcr1=1/sqrt(2*m)*sqrt(Ktot(1,1)-Ktot(2,2)/delta+sqrt((Ktot(1,1)+Ktot(2,2)/delta)^2-4*Ktot(1,2)^2/delta));
wcr2=1/sqrt(2*m)*sqrt(Ktot(1,1)-Ktot(2,2)/delta-sqrt((Ktot(1,1)+Ktot(2,2)/delta)^2-4*Ktot(1,2)^2/delta)); %complex

disp('Critical speed without gravity field [rad/s] = ')
disp(wcr1)
disp(' ')

figure
for j=1:size(Omega,2)
    wN(:,j)=roots([1 -Omega(j)*Jp/Jt -(Ktot(1,1)/m+Ktot(2,2)/Jt) +Omega(j)*Ktot(1,1)*Jp/(m*Jt) (Ktot(1,1)*Ktot(2,2)-Ktot(1,2)^2)/(m*Jt)]);
end
Omega_star=Omega/wcr1;
w_star=wN/wcr1;
plot(Omega_star,real(w_star),Omega_star,Omega_star,1,Omega_star,Omega_star,Jp/Jt*Omega_star)


    
   
