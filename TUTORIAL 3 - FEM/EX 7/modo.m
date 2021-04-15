function[fi,om_m]=modo(K,M,Nm,O)
%Nm e Ns numero gdl master e slaves

K11=K(1:Nm,1:Nm);
K12=K(1:Nm,Nm+1:end);
K21=K(Nm+1:end,1:Nm);
K22=K(Nm+1:end,Nm+1:end);

M11=M(1:Nm,1:Nm);

K_cond=K11-(K12*inv(K22))*K21;

%% solve
[fi_m,om_m]=eig(K_cond,M11);

fi_s=-(inv(K22)*K21)*fi_m;

fi=[fi_m; fi_s];

fi=riordinatore(fi,O,2);

