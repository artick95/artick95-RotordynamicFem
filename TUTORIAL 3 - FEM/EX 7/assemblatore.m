function[K]=assemblatore(plane)
%plane=1 if xz
%plane=2 if xy


%%matrices
%% DATI %%%
l1=0.3;
l2=0.123;
H=0.1; %altezza esterna beam
h=0.094; %altezza interna beam
alfa=deg2rad(20);
l=l1/cos(alfa);
theta=deg2rad(160);
D=0.035; %diametro esterno asta
s=0.001; %spessore asta
d=D-2*s;
E=7.2e10; %Pa;

%% Altre caratteristiche geometriche
At=H^2-h^2; %sezione beam
I=H^4/12-h^4/12;
Aa=pi*(D^2-d^2)/4;

%% Matrici
K_xz_12=E*[ At/l1    0            0    -At/l1       0           0;
    0       12*I/l1^3 6*I/l1^2     0     -12*I/l1^3 6*I/l1^2;
    0        6*I/l1^2  4*I/l1      0     -6*I/l1^2  2*I/l1;
    -At/l1      0        0      At/l1        0        0;
    0      -12*I/l1^3 -6*I/l1^2   0      12*I/l1^3 -6*I/l1^2;
    0        6*I/l1^2   2*I/l1    0      -6*I/l1^2    4*I/l1];


K_xz_23=E*[ At/l2    0            0    -At/l2       0           0;
    0       12*I/l2^3 6*I/l2^2     0     -12*I/l2^3 6*I/l2^2;
    0        6*I/l2^2  4*I/l2      0     -6*I/l2^2  2*I/l2;
    -At/l2      0        0      At/l2        0        0;
    0      -12*I/l2^3 -6*I/l2^2   0      12*I/l2^3 -6*I/l2^2;
    0        6*I/l2^2   2*I/l2    0      -6*I/l2^2    4*I/l2];


 
K_xy_12=E*[ 12*I/l1^3 6*I/l1^2          -12*I/l1^3 6*I/l1^2;
            6*I/l1^2  4*I/l1           -6*I/l1^2  2*I/l1;
          -12*I/l1^3 -6*I/l1^2         12*I/l1^3 -6*I/l1^2;
            6*I/l1^2   2*I/l1          -6*I/l1^2    4*I/l1];

K_xy_23=E*[12*I/l2^3  6*I/l2^2          -12*I/l2^3 6*I/l2^2;
           6*I/l2^2   4*I/l2            -6*I/l2^2  2*I/l2;
          -12*I/l2^3  -6*I/l2^2         12*I/l2^3 -6*I/l2^2;
           6*I/l2^2    2*I/l2          -6*I/l2^2    4*I/l2];        

        
cs=cos(theta);

sn=sin(theta);

K_xz_24=E*Aa/l*[cs^2   cs*sn  -cs^2 -cs*sn;
    cs*sn  sn^2   -cs*sn -sn^2;
    -cs^2  -cs*sn  cs^2 cs*sn;
    -cs*sn -sn^2   cs*sn sn^2];




%%maps
mappa_xz=[ 1 2 3 4 5 6 0 0 0 0 0;
    0 0 0 1 2 3 4 5 6 0 0;
    0 0 0 3 4 0 0 0 0 1 2];

mappa_xy=[1 2 3 4 0 0;
    0 0 1 2 3 4];

%%assemblying
if plane==1
    nel=3;
    Kbeam1=K_xz_12;
    Kbeam2=K_xz_23;
    Kbar=K_xz_24;
    mappa=mappa_xz;
else
    nel=2;
    Kbeam1=K_xy_12;
    Kbeam2=K_xy_23;
    mappa=mappa_xy;
end

Ngd=size(mappa,2);
K=zeros(Ngd);

for el=1:nel
    mapV=mappa(el,:);
    if el==1
        K_el=Kbeam1;
    end
    if el==2
        K_el=Kbeam2;
    end
    if el==3
        K_el=2*Kbar;
    end
    for ngd=1:Ngd
        if mapV(ngd)~=0
            for ind=ngd:Ngd
                if mapV(ind)~=0
                    K(ngd,ind)=K(ngd,ind)+K_el(mapV(ngd),mapV(ind));
                end
            end
        end
    end
end



