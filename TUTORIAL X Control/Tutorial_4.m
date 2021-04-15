%% 
% Dynamic Design of Machines
% Mechanical Engineering
% a.a. 2020-2021
% Tutorial 4  - Active control systems
% Ex. 1
%

%% Ex 1 - 1 dof

clear all
close all
clc

s               = tf('s');                                                  % Laplace variable  

m               = 1;                                                        % [kg]   - Mass
k               = 5e3;                                                      % [N/m]  - Stiffness of the spring k
k1              = 1e3;                                                      % [N/m]  - Stiffness of the spring k1
c               = 1;                                                        % [Ns/m] - Damping of the dashpot c
F               = 1;                                                        % [N]    - External force
Ksn             = 1;                                                        % [-]    - Sensor gain
Ccr             = 2*sqrt(k*m);

%c= Ccr;
KP              = ;                                                       % GR: 0.1; GS: 10
KD              = ;                                                      % GR: 0.1; GS: 0.5
KI              = ; 
Nd              = ;

% Gc_P          = KP;                                                       % Proportional controller
% Gc_PI         = KP*(1+1/(KI*s));                                          % Proportional Integrative controller
% Gc_PD         = KP*(1+(KD*s)/((KD/Nd)*s+1));                              % Proportional Derivative controller
% Gc_PID        = KP*(1+(1/(KI*s))+(KD*s)/((KD/Nd)*s+1));                   % Proportional Integrative Derivative controller
Gc              = KP + KD*s/(0.0001*s+1) + KI/s;                                         % [-]    - PID Control 
Gp              = 1/(m*s^2 + c*s + k);                                      % [m/N]  - Transfer function of the plant

figure, bode(Gp)
figure, step(Gp)
figure, pzmap(Gp)
pole(Gp)

GR              = 1/(m*s^2 + c*s + k + Gc*k1 + k1);                         % [m/N]  - Transfer function of the controlled system as regulator 
%GR             = minreal(GR,0.1);

figure, bode(Gp, GR)
figure, step(GR)
figure, pzmap(GR)
pole(GR)

GS              = (Gc*k1 + k1)/(m*s^2 + c*s + k + Gc*k1 + k1);              % [-]    - Transfer function of the controlled system as servomechanics
GS              = minreal(GS,0.1);

pole(GS)
figure, bode(Gp, GS)
figure, step(GS)
figure, pzmap(GS)