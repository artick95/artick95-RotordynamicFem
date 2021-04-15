clear all
close all
clc

%%%%%%%%%%%%%dire in che nodo si vuole lo stress%%%%%%%%%%%%%
nodo=1;
disp(['!!!!!!!!!!!!!!!!!' num2str(nodo) '!!!!!!!!!!!!!']);

g=9.81;
m=15;  %massa sospesa
E=7.2e10; %Pa

K_xz=assemblatore(1);
K_xy=assemblatore(2);

%% Matrici globali 
 
M_xz=zeros(11,11);
M_xz(7,7)=m;
M_xz(8,8)=m;

%manca M_xy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_xy=zeros(6,6);
M_xy(5,5)=m;
%%%%%%%%%%%%%%imporre i vincoli sulle matrici (eliminare alcune righe e
%%%%%%%%%%%%%%alcune colonne

M_xz= M_xz(4:9,4:9);
M_xy= M_xy(3:6,3:6);
K_xz= K_xz(4:9,4:9);
K_xy= K_xy(3:6,3:6);

M_xz= M_xz+(triu(M_xz,1)).';
M_xy= M_xy+(triu(M_xy,1)).';
K_xz= K_xz+(triu(K_xz,1)).';
K_xy= K_xy+(triu(K_xy,1)).';

%% riordino
Oxz=[4 1;
     5 2];

Oxy=[3 1];

M_xzr=riordinatore(M_xz,Oxz,1);
M_xyr=riordinatore(M_xy,Oxy,1);
K_xzr=riordinatore(K_xz,Oxz,1);
K_xyr=riordinatore(K_xy,Oxy,1);

%con la funzione "modo" estraggo autovalori e autovettori (ordinati non
%secondo slaves e master, ma come erano ordinati all'inizio)
Nmxz=2;       %numero master in xz
Nmxy=1;      %numero master in xy


[fi_xz,auto_xz]=modo(K_xzr,M_xzr,Nmxz,Oxz);
[fi_xy,auto_xy]=modo(K_xyr,M_xyr,Nmxy,Oxy);

om_xz=diag(auto_xz.^0.5);
om_xy=diag(auto_xy.^0.5);

%% risposta forzata armonica

%modali
M_xz_mod=diag(fi_xz.'*M_xz*fi_xz);
M_xy_mod=diag(fi_xy.'*M_xy*fi_xy);

%normalizzazione
fi_xz(:,1)=fi_xz(:,1)/sqrt(M_xz_mod(1));
fi_xz(:,2)=fi_xz(:,2)/sqrt(M_xz_mod(2));
fi_xy(:,1)=fi_xy(:,1)/sqrt(M_xy_mod(1));

delta_x=[1;0;0;1;0;0];

delta_z=[0;1;0;0;1;0];

delta_y=[1;0;1;0];

eta_xz=@(lamb,Xa0,Za0) (-lamb.^2*eye(2)+auto_xz)\(lamb.^2*fi_xz.'*M_xz*(delta_x*Xa0+delta_z*Za0));
eta_01=@(lamb,Ya0) lamb.^2*(fi_xy(:,1)'*M_xy*delta_y*Ya0)./(om_xy(1)^2-lamb.^2);

df=0.1; %risoluzione in frequenza
dim=(50-5)/df+1;
ind=1;
Qxz=zeros(6,dim); %sulle colonne si spazza la frequenza, sulle righe si hanno i DOF
Qxy=zeros(4,dim);

%out in m
for fr=5:df:8.5-df
    f=fr*2*pi;
    Xa0=0.010; %in m
    Za0=Xa0;
    Ya0=Xa0;
    n_01=eta_01(f,Ya0);
    eta_p=eta_xz(f,Xa0,Za0);
    n_02=eta_p(1);
    n_03=eta_p(2);
    Qxz(:,ind)=fi_xz*[n_02; n_03];
    Qxy(:,ind)=fi_xy*n_01;
    ind=ind+1;
end

for fr=8.5:df:35-df
    f=fr*2*pi;
    Xa0=3*g/f^2;
    Za0=Xa0;
    Ya0=Xa0;
    n_01=eta_01(f,Ya0);
    eta_p=eta_xz(f,Xa0,Za0);
    n_02=eta_p(1);
    n_03=eta_p(2);
    Qxz(:,ind)=fi_xz*[n_02; n_03];
    Qxy(:,ind)=fi_xy*n_01;
    ind=ind+1;
end

for fr=35:df:50
    f=fr*2*pi;
    Xa0=g/f^2;
    Za0=Xa0;
    Ya0=Xa0;
    n_01=eta_01(f,Ya0);
    eta_p=eta_xz(f,Xa0,Za0);
    n_02=eta_p(1);
    n_03=eta_p(2);
    Qxz(:,ind)=fi_xz*[n_02; n_03];
    Qxy(:,ind)=fi_xy*n_01;
    ind=ind+1;
end

[U3x_max,induxmax]=max(Qxz(4,:));
[U3y_max,induymax]=max(Qxy(3,:));
[U3z_max,induzmax]=max(Qxz(5,:));


for ind=1:dim
    [MyI(ind),MyII(ind),My3II(ind),MzI(ind),MzII(ind),Mz3II(ind),sigma(ind)]=forze(Qxz(:,ind),Qxy(:,ind),nodo);
end

%sollecitazioni al nodo 1 o 2 struttura
if nodo==2
    figure;
    plot(5:df:50,MyI,5:df:50,MyII);
    title('My, Harmonic');
    legend('ElemI','ElemII');
    xlabel('frequency [Hz]','Interpreter','LaTex');
    ylabel('My [Nm]','Interpreter','LaTex');
    figure;
    plot(5:df:50,MzI,5:df:50,MzII)
    title('Mz, Harmonic');
    legend('ElemI','ElemII');
    xlabel('frequency [Hz]','Interpreter','LaTex');
    ylabel('Mz [Nm]','Interpreter','LaTex');
elseif nodo==1
    figure;
    plot(5:df:50,MyI);
    title('My, Harmonic');
    legend('ElemI node 1');
    xlabel('frequency [Hz]','Interpreter','LaTex');
    ylabel('My [Nm]','Interpreter','LaTex');
    figure;
    plot(5:df:50,MzI)
    title('Mz, Harmonic');
    legend('ElemI node 1');
    xlabel('frequency [Hz]','Interpreter','LaTex');
    ylabel('Mz [Nm]','Interpreter','LaTex');    
end

%stress
%solo su elemento 1 al nodo 1 oppure 2
figure; 
plot(5:df:50,sigma*1e-6);
title('Stress, Harmonic');
xlabel('frequency [Hz]','Interpreter','LaTex');
ylabel('$\sigma$ [MPa]','Interpreter','LaTex');

%% shock
%parametri in simulink
lamb1=om_xy(1);
Ry=fi_xy.'*M_xy*delta_y;

%durata simulazione
tfin=0.2;
sim Shock_sim
q_shock=q_shock';
%in out da modello simulink vettore q_shock, prima riga è il vettore dei
%tempi; sulle altre righe ci sono i DOF, lungo le colonne si spazzano i
%tempi

dim=size(q_shock,2);
for ind=1:dim
    [MzI_shock(ind),MzII_shock(ind),Mz3II_shock(ind),sigma_shock(ind)]=forze_shock(q_shock(2:end,ind),nodo);
end

U3y_shock_max=max(q_shock(4,:));
if nodo==2
figure; 
plot(q_shock(1,:),MzI_shock,q_shock(1,:),MzII_shock)
title('Mz, shock');
legend('ElemI','ElemII');
xlabel('frequency [Hz]','Interpreter','LaTex');
ylabel('Mz [Nm]','Interpreter','LaTex');
elseif nodo==1
    figure;
    plot(q_shock(1,:),MzI_shock)
    title('Mz, shock');
    legend('ElemI');
    xlabel('frequency [Hz]','Interpreter','LaTex');
    ylabel('Mz [Nm]','Interpreter','LaTex');
end
%stress
%solo su elemento 1 al nodo 2
figure; 
plot(q_shock(1,:),sigma_shock*1e-6);
title('Stress, shock');
xlabel('Time [s]','Interpreter','LaTex');
ylabel('$\sigma$ [MPa]','Interpreter','LaTex');

%% verifiche
Rm=328e6; %Pa
Rp02=216e6;%Pa
sd_1=115e6; %Pa, limite di fatica

%HARMONIC
Umax=max([U3x_max,U3y_max,U3z_max]);
disp('HARMONIC');
disp(['max disp in mm ' num2str(Umax*1e3)]);

if max(sigma)<=Rm/1.575
    disp(['Rottura verificata, sigma [Mpa]=' num2str(max(sigma/1e6))]);
else
    disp(['Rottura non verificata, sigma  [Mpa]=' num2str(max(sigma/1e6))]);
end

if max(sigma)<=Rp02/1.55
    disp(['Snervamento verificato, sigma  [Mpa]=' num2str(max(sigma/1e6))]);
else
    disp(['Snervamento non verificato, sigma [Mpa]=' num2str(max(sigma/1e6))]);
end

if max(sigma)<=sd_1
    disp(['Fatica verificata, sigma  [Mpa]=' num2str(max(sigma/1e6))]);
else
    disp(['Fatica non verificata,, sigma [Mpa]=' num2str(max(sigma/1e6))]);
end

%IMPULSE
disp('SHOCK');
disp(['max disp in mm ' num2str(U3y_shock_max*1e3)]);

if max(sigma_shock)<=Rm
    disp(['Rottura verificata, sigma [Mpa]=' num2str(max(sigma_shock/1e6))]);
else
    disp(['Rottura non verificata, sigma [Mpa]=' num2str(max(sigma_shock/1e6))]);
end












