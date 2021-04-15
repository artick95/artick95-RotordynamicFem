clear all
clc

E=2e11;
l=.125;
d=.025;
m1=30;
m2=50;
L=5*l;

I=pi*d^4/64;
a1=l;
b1=4*l;
x1=(a1*b1)^2/(3*E*I*L);
a2=2*l;
b2=3*l;
x2=(a2*b2)^2/(3*E*I*L);
wcr=sqrt(1/(x1*m1+x2*m2));

disp('First bending Critical speed [rad/s] =')
disp(wcr)