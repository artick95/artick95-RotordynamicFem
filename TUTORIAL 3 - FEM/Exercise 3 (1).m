close all
clear
clc

% unit [mks]

P = 10000;
LL = 1;
E = 210e9;

% element I
L1 = LL;
h1 = 30e-3;
b1 = 15e-3;
A1 = h1*b1;
I1 = b1*h1^3/12;
alpha1 = pi/4;

% element II
L2 = LL;
h2 = 20e-3;
b2 = 10e-3;
A2 = h2*b2;
I2 = b2*h2^3/12;
alpha2 = pi-pi/4;


%% STRUCTURE 3A
Rz =@(t) [cos(t) sin(t); -sin(t) cos(t)];
Rz = @(t) [Rz(t) zeros(2); zeros(2) Rz(t)];

A = A1;
L = L1;
Kt1 = E*A/L*[1  0   -1  0;
            0   0   0   0;
            -1  0   1   0;
            0   0   0   0];
	
Kt1g = (Rz(alpha1)')*Kt1*Rz(alpha1);

A = A2;
L = L2;
Kt2 = E*A/L*[1  0   -1  0;
            0   0   0   0;
            -1  0   1   0;
            0   0   0   0];
	
Kt2g = (Rz(alpha2)')*Kt2*Rz(alpha2);

% assembly
Kg = zeros(6);
Kg(1:4,1:4) = Kt1g;
Kg(3:6,3:6) = Kg(3:6,3:6)+Kt2g;

% boundary conditions
Kred = Kg([3 4],[3 4]);
Fext = [0; -P];

usol = Kred\Fext;

uvec = [0; 0; usol; 0; 0];
Fvec = Kg*uvec;


%% STRUCTURE 3B
Rz =@(t) [cos(t) sin(t) 0 ; -sin(t) cos(t) 0; 0 0 1];

Rz = @(t) [Rz(t) zeros(3); zeros(3) Rz(t)];

A = A1;
L = L1;
I = I1;
Kb1 = E*[A/L    0           0           -A/L        0           0;
        0       12*I/L^3    6*I/L^2     0           -12*I/L^3   6*I/L^2;
        0       6*I/L^2     4*I/L       0           -6*I/L^2    2*I/L;
        -A/L     0           0           A/L        0           0;
        0       -12*I/L^3   -6*I/L^2    0           12*I/L^3    -6*I/L^2;
        0       6*I/L^2     2*I/L       0           -6*I/L^2    4*I/L;
        ];
	
Kb1g = (Rz(alpha1)')*Kb1*Rz(alpha1);

A = A2;
L = L2;
I = I2;
Kb2 = E*[A/L    0           0           -A/L        0           0;
        0       12*I/L^3    6*I/L^2     0           -12*I/L^3   6*I/L^2;
        0       6*I/L^2     4*I/L       0           -6*I/L^2    2*I/L;
        -A/L     0           0           A/L        0           0;
        0       -12*I/L^3   -6*I/L^2    0           12*I/L^3    -6*I/L^2;
        0       6*I/L^2     2*I/L       0           -6*I/L^2    4*I/L;
        ];
	
Kb2g = (Rz(alpha2)')*Kb2*Rz(alpha2);

% assembly
Kg = zeros(10);
Kg(1:6,1:6) = Kb1g;
Kg([4 5],[4 5]) = Kg([4 5],[4 5])+Kb2g([1 2],[1 2]);
Kg([4 5],7:10) = Kg([4 5],7:10)+Kb2g([1 2],3:6);
Kg(7:10,[4 5]) = Kg(7:10,[4 5])+Kb2g(3:6,[1 2]);
Kg(7:10,7:10) = Kg(7:10,7:10)+Kb2g(3:6,3:6);

% boundary conditions
Kred = Kg(4:7,4:7);
Fext = [0; -P; 0; 0];

usol = Kred\Fext;

uvec = [0; 0; 0; usol; 0; 0; 0];
Fvec = Kg*uvec;


%% STRUCTURE 3C
Rz =@(t) [cos(t) sin(t) 0 ; -sin(t) cos(t) 0; 0 0 1];

Rz = @(t) [Rz(t) zeros(3); zeros(3) Rz(t)];

A = A1;
L = L1;
I = I1;
Kb1 = E*[A/L    0           0           -A/L        0           0;
        0       12*I/L^3    6*I/L^2     0           -12*I/L^3   6*I/L^2;
        0       6*I/L^2     4*I/L       0           -6*I/L^2    2*I/L;
        -A/L     0           0           A/L        0           0;
        0       -12*I/L^3   -6*I/L^2    0           12*I/L^3    -6*I/L^2;
        0       6*I/L^2     2*I/L       0           -6*I/L^2    4*I/L;
        ];
	
Kb1g = (Rz(alpha1)')*Kb1*Rz(alpha1);

A = A2;
L = L2;
I = I2;
Kb2 = E*[A/L    0           0           -A/L        0           0;
        0       12*I/L^3    6*I/L^2     0           -12*I/L^3   6*I/L^2;
        0       6*I/L^2     4*I/L       0           -6*I/L^2    2*I/L;
        -A/L     0           0           A/L        0           0;
        0       -12*I/L^3   -6*I/L^2    0           12*I/L^3    -6*I/L^2;
        0       6*I/L^2     2*I/L       0           -6*I/L^2    4*I/L;
        ];
	
Kb2g = (Rz(alpha2)')*Kb2*Rz(alpha2);

% assembly
Kg = zeros(9);
Kg(1:6,1:6) = Kb1g;
Kg(4:9,4:9) = Kg(4:9,4:9)+Kb2g;

% boundary conditions
Kred = Kg(4:6,4:6);
Fext = [0; -P; 0];

usol = Kred\Fext;

uvec = [0; 0; 0; usol; 0; 0; 0];
Fvec = Kg*uvec;