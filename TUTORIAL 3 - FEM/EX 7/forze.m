function[MyI,MyII,My3II,MzI,MzII,Mz3II,sigma]=forze(Qxz,Qxy,nodo)

qIxz=[0; 0; 0; Qxz(1:3)];
qIIxz=Qxz;
qIxy=[0; 0; Qxy([1;2])];
qIIxy=Qxy;

%% matrices
%DATI
l1=0.3;
l2=0.123;
H=0.1; %altezza esterna beam
h=0.094; %altezza interna beam

E=7.2e10; %Pa;

%Altre caratteristiche geometriche
At=H^2-h^2; %sezione beam
I=H^4/12-h^4/12;

%Matrici
KIxz=E*[ At/l1    0            0    -At/l1       0           0;
    0       12*I/l1^3 6*I/l1^2     0     -12*I/l1^3 6*I/l1^2;
    0        6*I/l1^2  4*I/l1      0     -6*I/l1^2  2*I/l1;
    -At/l1      0        0      At/l1        0        0;
    0      -12*I/l1^3 -6*I/l1^2   0      12*I/l1^3 -6*I/l1^2;
    0        6*I/l1^2   2*I/l1    0      -6*I/l1^2    4*I/l1];


KIIxz=E*[ At/l2    0            0    -At/l2       0           0;
    0       12*I/l2^3 6*I/l2^2     0     -12*I/l2^3 6*I/l2^2;
    0        6*I/l2^2  4*I/l2      0     -6*I/l2^2  2*I/l2;
    -At/l2      0        0      At/l2        0        0;
    0      -12*I/l2^3 -6*I/l2^2   0      12*I/l2^3 -6*I/l2^2;
    0        6*I/l2^2   2*I/l2    0      -6*I/l2^2    4*I/l2];


 
KIxy=E*[ 12*I/l1^3 6*I/l1^2          -12*I/l1^3 6*I/l1^2;
            6*I/l1^2  4*I/l1           -6*I/l1^2  2*I/l1;
          -12*I/l1^3 -6*I/l1^2         12*I/l1^3 -6*I/l1^2;
            6*I/l1^2   2*I/l1          -6*I/l1^2    4*I/l1];

KIIxy=E*[12*I/l2^3  6*I/l2^2          -12*I/l2^3 6*I/l2^2;
           6*I/l2^2   4*I/l2            -6*I/l2^2  2*I/l2;
          -12*I/l2^3  -6*I/l2^2         12*I/l2^3 -6*I/l2^2;
           6*I/l2^2    2*I/l2          -6*I/l2^2    4*I/l2];
       
%% 

FIxz=KIxz*qIxz;
FIIxz=KIIxz*qIIxz;
FIxy=KIxy*qIxy;
FIIxy=KIIxy*qIIxy;

%% sollecitazioni
if nodo==2
    FxI=FIxz(4,:);
    MyI=FIxz(end,:);
    MyII=FIIxz(3,:);
    My3II=FIIxz(end,:);
    MzI=FIxy(end,:);
    MzII=FIIxy(2,:);
    Mz3II=FIIxy(end,:);
    
elseif nodo==1
    FxI=FIxz(1,:);
    MyI=FIxz(3,:);
    MyII=0;
    My3II=0;
    MzI=FIxy(2,:);
    MzII=0;
    Mz3II=0;
    
end
    

%%stress
sigma=abs(MyI)/I*H/2+abs(MzI)/I*H/2+abs(FxI/At);