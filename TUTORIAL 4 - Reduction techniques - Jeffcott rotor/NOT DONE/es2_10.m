clear all
close all
clc

E=2e11;
v=.3;
L=.1;
d=.03;
m=2;
Jt=.005;
Jp=.010;

G=E/(2*(1+v));
A=pi*d^2/4;
I=pi*d^4/64;
X=(7+6*v)/(6*(1+v));
phi=12*E*I*X/(G*A*L^2);

wcr_s=input(' Guess a critical speed considering shear stress : ');
disp(' ')
wcr=input(' Guess a critical speed NOT considering shear stress (it should be higher): ');
disp(' ')

Tf_s=[1   L   -L^3/(6*E*I)+L*X/(G*A)    L^2/(2*E*I)
      0   1   -L^2/(2*E*I)              L/(E*I)
      0   0   1                         0
      0   0   -L                        1];

Tf=[1   L   -L^3/(6*E*I)            L^2/(2*E*I)
    0   1   -L^2/(2*E*I)            L/(E*I)
    0   0   1                       0
    0   0   -L                      1];

Tn_s=[1           0                   0   0
      0           1                   0   0
      -wcr_s^2*m  0                   1   0
      0           -wcr_s^2*(Jt-Jp)    0   1];

Tn=[1           0                   0   0
    0           1                   0   0
    -wcr^2*m    0                   1   0
    0           -wcr^2*(Jt-Jp)      0   1];

% Consider shear                
T_s=Tn_s*Tf_s*eye(4);                
err_s=T_s(3,3)*T_s(4,4)-T_s(4,3)*T_s(3,4); 
disp('Critical speed considering shear deformation [rad/s] =')
disp(wcr_s)
disp(' ')
disp('Error =')
disp(err_s)
disp(' ')

% Not considering shear (X=0 hence phi=0) 
T=Tn*Tf*eye(4);                
err=T(3,3)*T(4,4)-T(4,3)*T(3,4); 
disp('Critical speed not considering shear deformation [rad/s] =')
disp(wcr)
disp(' ')
disp(' ')
disp('Error =')
disp(err)
disp(' ')


% Using REST PLOT

% Consider shear  
hold on
fplot(@(w) (w^2*m*(L^3/(6*E*I)-L*X/(G*A))+1)*(-w^2*(Jt-Jp)*L/(E*I)+1)-(-w^2*m*L^2/(2*E*I))*(w^2*(Jt-Jp)*L^2/(2*E*I)-L),[0 10000])
plot(linspace(0,10000,10001),0)
plot(4.294e3,linspace(-1,1,2001))

% NOT Consider shear 
figure
hold on
fplot(@(w) (w^2*m*(L^3/(6*E*I))+1)*(-w^2*(Jt-Jp)*L/(E*I)+1)-(-w^2*m*L^2/(2*E*I))*(w^2*(Jt-Jp)*L^2/(2*E*I)-L),[0 10000])
plot(linspace(0,10000,10001),0)
plot(4.552e3,linspace(-1,1,2001))
