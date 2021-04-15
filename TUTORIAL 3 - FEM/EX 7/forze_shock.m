function[MzI,MzII,Mz3II,sigma]=forze_shock(Qxy,nodo)

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
KIxy=E*[ 12*I/l1^3 6*I/l1^2          -12*I/l1^3 6*I/l1^2;
            6*I/l1^2  4*I/l1           -6*I/l1^2  2*I/l1;
          -12*I/l1^3 -6*I/l1^2         12*I/l1^3 -6*I/l1^2;
            6*I/l1^2   2*I/l1          -6*I/l1^2    4*I/l1];

KIIxy=E*[12*I/l2^3  6*I/l2^2          -12*I/l2^3 6*I/l2^2;
           6*I/l2^2   4*I/l2            -6*I/l2^2  2*I/l2;
          -12*I/l2^3  -6*I/l2^2         12*I/l2^3 -6*I/l2^2;
           6*I/l2^2    2*I/l2          -6*I/l2^2    4*I/l2];
       
%% 
FIxy=KIxy*qIxy;
FIIxy=KIIxy*qIIxy;

%% sollecitazioni
if nodo==2
    MzI=FIxy(end,:);
    MzII=FIIxy(2,:);
    Mz3II=FIIxy(end,:);
    
elseif nodo==1
    MzI=FIxy(2,:);
    MzII=0;
    Mz3II=0;
end

%%stress
sigma=(MzI)/I*H/2;
