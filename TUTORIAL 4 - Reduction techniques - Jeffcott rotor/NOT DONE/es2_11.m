clear all
clc

d=.04;
l=.6;
gamma=30;
E=2e11;
m=5;

r=d/2;

% Consider only half circle, subtract the later circular sector and add the
% triangles
% Ie_c=pi*r^4/8
% Ie_cs=r^4/8*(gamma*pi/90-sin(gamma*pi/90))
% Ie_t=r^4/12*cosd(gamma)*(sind(gamma))^3
% Ie=2*(Ie_c-Ie_cs+2*Ie_t)
Ie=2*(pi*r^4/8-r^4/8*(gamma*pi/90-sin(gamma*pi/90))+r^4/6*cosd(gamma)*(sind(gamma))^3);

% Consider only one circular sector and only one triangle:
% In_cs=r^4/8*((180-2*gamma)*pi/180-sin((180-2*gamma)*pi/180))
% In_t=r^4/4*sind(gamma)*(cosd(gamma))^3;
% In=2*In_cs+4*In_t
In=r^4/4*((180-2*gamma)*pi/180-sin((180-2*gamma)*pi/180))+r^4*sind(gamma)*(cosd(gamma))^3;

ke=48*E*Ie/l^3;
kn=48*E*In/l^3;

wcr1=sqrt(ke/m);
wcr2=sqrt(kn/m);

disp('Range of instability [rad/s] =')
disp(' ')
disp(wcr2) 
disp(wcr1)