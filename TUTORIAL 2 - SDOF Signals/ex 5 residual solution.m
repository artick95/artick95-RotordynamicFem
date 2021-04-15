
clear
s=tf('s');

T0=5/1000;
F0=10;
m=2;
k = 100;
c=0;

H=1/(m*s^2+c*s+k);
F=F0/s^2;
X=H*F;
zpk(minreal(X))

[Num,Den]=tfdata(X,'v');
[r,p]=residue(Num,Den)