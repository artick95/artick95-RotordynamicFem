stop = 1

while (stop != 0):
    print(
        'how can I help?\n [1]SDOF\n [2]Responses\n [3]Intro+Pr\n [4]Control Systems\n [5]MDOF\n [6]Signals\n [7]Laplace\n [8]FEM\n [9]Rotordynamics\n [10]Reduction Tecniques\n [11]Torsional Vibrations\n [12]FRFz\n [13]MathStuff')

    formulaSheets = ['You fucked', 'SDOF&MDOF','Responses','Intro','Control Systems','MDOF','Signals','Laplace','FEM','Rotordynamics','Reduction Tecniques','Torsional Vibrations','FRFz','MathStuff']

    formulaSheet_chosen = int(input())

    print(formulaSheets[formulaSheet_chosen])



    if (formulaSheet_chosen == 1):
        print(' **SDOF**\n Ka=A*E/l\n Kb=EI/l^3\n Kt=G*Ip/l\n Kshear=12G*A/(l*xsi)\n Kbcantilever=3EI/l^3\n G=7.7GPa\n E=210GPa\n taumax=Mt*d/2/Ip\n **Inertia Proprieties**\n Wt=2Wf\n Wf=pi*d^3/32\n Ip=2Ixx\n Ixx=pi*d^4/64\n Ix=b*h^3/12\n Mf=-EIy,,\n Mt = G*I*teta.\n x+y+z+x+\n gammax = xsi*T/(G*A)\n phi=EI*l*xsi/(l^3*pi*G*A)\n **Hysteretic damping** \n ceq = eta*k/omega\n eta = 2*zeta\n phi=2*eta=ellipse/triangle0AB\n c=c-+ic..\n E=E.+iE..\n k=k.+k..\n eta=phi=atan(phi)\n k,=k\n Viscoud damping capacity PHI=Ed/Es\n k*=k(1+ieta)\n Ed=pi*c*Omega*x0^2=pisigmazeroepsilon0sinphi \n **Avarage Eta Values**\n 0.00001%0.001 Al\n 0.01%0.06 steel\n 0.01%1,3 Rubber\n **Avarage Zeta Values**\n 0.001 material \n 0.2%0.3 elective active system\n 0.03 cast iron')
        extra = 0
        print('load more?[1-y][0-n]')
        extra = int(input())
        if extra == 1:
            print(' **Logaritimic DecMet**\n delta=ln(xi/x(i+1))=2*pi*zita=2*pi*zita/sqrt(1-zita^2)\n zita=delta/(2*pi)=c/cd\n  **Maxwerl-Wircht**\n F=-kx-int(C(t,tau,x.(tau)dtau,-inf,t)\n C(t,tau)=sum(c_imu_ie^(-mu_i(t-tau)))\n k(s)=ci*mui*s/(mui+s)')

    if (formulaSheet_chosen == 2):
        print(' **Responses**\n **Genericnonstaitonaryinput**\n x=K*e^^(-zetawnt)*e^(iwnsqrt(1-zita^2)*t)+H(omega)*f0/k*e^(iwt);K=x0\n **Impulse**\n x=f0/(m*wn)*h\n h=1/sqrt(1-zeta^2))*e^(-zetawn*t)*sin(wn*sqrt(1-zeta^2)*t)\n h=wn*t*e^(-wnt)\n h=1/(2*sqrt(1-zita^2)){-e^(-(zeta+sqrt(1-zeta^2))*wnt)+e^(-(zeta-sqrt(1-zeta^2))*wnt)}\n **Step**\n x=f0/k*g\n g=1-e^-(zetawnt)*[cos(wn*sqrt(1-zeta^2)*t)+zeta/sqrt(1-zeta^2)*sin(wn*sqrt(1-zeta^2)*t)\n g=1-(1-wn*t)*e^(wn*t)\n g=1-0.5*{-e^(-(zeta+sqrt(1-zeta^2))*wn*t)+e^(-(zeta-sqrt(1-zeta^2))*wn*t)}\n **Duhamel** \n x=1/(m*wn)*int(F(tau)*h(t-tau)dtau)\n x=1/(m*wn)*int(F(t)*h(tau)dtau) \n x=Asin(sqrt(1-zeta^2)*wn*t)-B*cos(sqrt(1-zeta^2)wn*t)\n A=1/(wn*m*sqrt(1-zeta^2))*e^(zeta*wn*t)int(F(tau)*e^(zetawntau)*cos(sqrt(1-zeta^2)*wn*tau)dtau\n B = A with sine  ')

    if (formulaSheet_chosen == 3):
        print(' **Intro+Pr**\n exposure limit *2\n reduced confort boundary/3.15\n 4-8Hz human more affected\n <1Hz Motion sickness\n 18-20 Hz audible intervall\n static = self weight of a building\n quasi-static = w<<wn (resonance).The loads are applied slowly and remain for a long period of time at more or less constant values. e.g. "Centrifugal force"\n dynamic = w=wn\n non linear system do not allow any transform\n **Project1**\n hardtip-deltaT smaller\n wcutoff=1/deltaT\n kdeltax=sensitivit*V\n ux0..=wn^2*sensitivity/k\n sensitivity:V/s\n C=phit^-1*Cmod*phi^-1\n zitamod=Cmod/Ccrit\n Copt=sqrt(kms/2)\n c_emteta=Kt^2/(Ra+Rext_variable)\n c_em_lin=c_em_teta*eta^2\n Sx(ni)=Gr*ni^(-2)\n m^2/(cycle/m)\n m*cycle\n (m/cycle)^2\n freq=speed*ni(space)')
        extra = 0
        print('load more?[1-y][0-n]')
        extra = int(input())
        if extra == 1:
            print(' **Project2**\n kdeltax=V*Sensitivity \n  **Project 3 Stuff**\n Kp=0 Td=0 unstable 2 poles simmetric to y\n  Td=0 KpKiKsn-Kx>0 Marginally stable 2 simmetric on x axis\n Td>0 KpKiKsn-Kx>0 Stable all to left 1 in origin \n **PID**\n uc= Kpey+Kd*dey/dt*Kiint(eydt) ')
    if (formulaSheet_chosen == 4):
        print(' **ControlSystems**\n Adaptive: actuator\n sensory: sensor\n controlled: adaptive+sensory\n active: controlled with reference input\n intelligent: bird\n Gyr=GcGaGp/(1+L) (+if-in feedback)\n L=G*H\n TcKpTs positive defined up stability ')

    if (formulaSheet_chosen == 5):
        print(' **MDOF**\n **Lagrange Equation**\n diff(diff(T,x.),t)-diff(T,x)+diff(UA,x)=Qi\n Mx..+(C+G)x.+(K+H)x=f\n M—>conservative\n C(symmetric linear dissipation)—>non conservative but if C is positive defined stabilizing\n G—>does not affect stability( circulatory/gyroscopic and coriolis system skew-symmetric)\n K—>does not affect stability(simmetric)\n H—>if H not positive defined can cause instability(gyroscopic system skew-symmetric)\n C and K —> symmetric\n Positive define:xT*A*x>0\n M not diagonal inertial coupling\n K not diagonal elastic coupling\n F=1/2x.TCx.+x.THx\n H screw symmetric\n G=0 and H=0 —> natural system\n G=H=C=0 —> MK system\n H=0 —> gyroscopic system skew-symmetric\n  G=0—> circulatory system skew-symmetric\n z={v;x}\n z.=Az+Bu\n y=Cz+Du\n F=0.5x.TCx.+x.THx\n if det(M)=0 singular ')

    if (formulaSheet_chosen == 6):
        print(' **Signals**\n white noise has all frequencies but yrms=inf(bad)\n A state vector identity of the system \n b input gain matrix forced +homogeneous behavior\n C output gain matrix\n D direct link matrix\n z. memory\n z\n A not depends on input only on MCGKH\n yrms=sqrt(int(S(w),-int,inf))\n MIL 20g 11ms\n Sout=Sin*|H|^2\n stationary: you can start the analysis at any time\n ergodic: there is a subset in the hystory that is rappresentative\n ymed=1/T*int(ydt)\n ymed2=1/T*int(y^2dt)\n yrms=sqrt(ymed2)\n PSD: yrms=sqrt(int(S(w)*dw)) -int +inf\n variance=sqrt(mean(y-ymed)^2)  ')

    if (formulaSheet_chosen == 7):
        print(' **Laplace**\n L(deltadirac)=1\n L(F0(ut-ut-5))=F0/s-F0/s*e^-T0s \n L(e^(lambda*t))=1/(s-lambda)\n L_1=e^(lambda*t)\n L(1)=1/s\n L_1=1\n L(t^n)=n!/s^(n+1)\n L_1(1/s^n)=t^(n-1)/(n-1)!\n L(cosbetat)=s/(s^2+beta^2)\n L_1(s/(s^2+beta^2))=cosbetat\n L(sinbetat)=beta/(s^2+beta^2)\n L_1(1/(s^2+beta^2))=1/beta*sin(beta*t)\n L(x.)=sL(x)-x0\n L(x..)=s^2*L(x)-s*x0-x.0\n *e^(-T0*s):shift to the right in laplace world\n L(f)=int(f*e^-st)dt) 0 to +inf\n L_1(F(s+alfa))=L_1(F(s)*e^-alfat)\n L_1(1-e^-5t)=heaviside(t)-heaviside(t-5)')

    if (formulaSheet_chosen == 8):
        print(' **FEM**\n K=int(AbTEb,z)=int(BtEBdV)\n m=int(rhoANtNdz)=int(rhoNtNdV\n M=rhoAl/6[2 1,1 2]\n epsilon=B*q(t)\n u(x,y,z,t)=N(x,y,y)q(t)\n sigma=EBq\n U=1/2int(epsilon^T*sigmadV)=1/0q^Tint(B^T*EBdV)q\n K=int(B^T*EBdV)\n T=q.^Tint(rhoN^TNdV)q.=rho*dV*A*0.5*q.^2\n M=int(rhoN^TNdV)\n **Kirchoff**\n phiy=-duz/dx\n phix=duz/dy\n 3dof*4nodi=12dof*element\n **N Req**\n -continuousC1,rigid mot const strain, deflection sim neiborgh\n sigma=sigma0+E(epsilon-epsilon0)\n **EulerBeam**\n phiy=dux/dz\n **Timoshenko**\n dux/dz=gammax+phiy\n  Fx=GA/xsi*gammax\n Nyz=Nxz cambiato di segno 21 12 23 14\n 1st formulation inconsistency:\n dphy/dz=const->M=const gamma 1st order\n 2nd formulation\n tau = const\n gammax=((ux2-ux1)/l-(phiy1+phiy2/2)phi/(phi+1)\n phi=(xsi/(GA)/(L^2/(12eiY)\n phi up N down \n phi=0 slender beam  ')

    if (formulaSheet_chosen == 9):
        print(' **Rotordynamics**\n mz..+(cr+cn)z.+(k-icrteta.)z=mepsilon(teta.^2-iteta..)e^(iteta)\n z0=epsilon*r^2/sqrt((1-r^2)^2+4zita_n^2r^2)\n z0=epsilon*Omega/(wcr^2-Omega^2)\n zitar=cr/(2sqrt(km))\n wth<sqrt(k/m)*(1+cn/cr)\n stability Condition: wI<0 & wR<0\n cn always stabilizing\n cr instabilizing for FWD w>wth\n **Composite Pendulum**\n lteta..+Omega^2*rteta=0\n lphi..+Omega^2(r+l)phi\n w=Omegasqrt(1+r/l)\n w2=Omegasqrt(r/l)(use this for examstipe)\n cincronicity Omega=teta.\n **Structural Damping**\n z0=epsilon*r^2/sqrt((1-r^2)+etan^2)\n phi=arctg(-etan/(1-r^2))\n cn>etar*kr*sqrt(m/k)=etar*kr/Omegacr shaft\n Kshaft=48EI/l^3\n delta=Jp/Jt\n disk:Jp>Jt 1 int\n beam:Jp<Jt 2 int\n Sphere:Jp=Jt @INF\n zetan=cn/cr\n zeta=ctot/ccritical\n zetar=cr/ccritical')
        extra=0
        print('load more?[1-y][0-n]')
        extra=int(input())
        if extra==1:
            print('[m 0,0 Jt]{z.. phi..}+[K11 K12,K12 K22]{z phi} - i[0 0,0 Jp]Omega{z. phi.}=Omega^2{m*epsilon*e^^ialfa) xsi*(Jt-Jp)}*e^(i*Omega*t)\n above Force\n below Moment\nWcr=1/sqrt(2m)*sqrt(K11-K22/delta+-sqrt((K11+K22/delta)^2-4K12^/delta))\n delta=(Jp-Jt)/m\n*Dipendenza a*\n w=OmegaJp/Jt \n Mg=-i*Jp*Omega*phi.+Jtphi..\n **NonIsJefcot**\n mx..+kx*x=mepsilon*omega^2*cos(omegat)\n my..+ky*y=mepsilon*omega^2*sin(omegat)\n  km=0.5(kx+ky)\n kd=0.5(kx-ky)\n mz..+km*z+kdzconj=mepsilonOmega^2e^(iOmegat)\n **Hydrodinamic Bearing**\n Fx,=0\n Fy,=12pi*mu*l*(Rj/c)^3*(Omega/2-omega)*epsilon\n **OilWhirlWhip**\n w<Omega/2 unstable(real0.45-0.48Omega)\n  **Balancing**\n Balancing grade=Omega*epsilon \n xsi(JP-jT)=mepsilond/2\n Rigid: Omega<Omegacr static balance\n Flex Field balancing  ')
            extra = 0
            print('load more?[1-y][0-n]')
            extra = int(input())
            if extra == 1:
                print('**4DOF**\n Mx=Jtphi..x,+JpOmegaphi.y+xsiOmega^2(Jt-Jp)sin(Omegat)\n My=Jtphi..y,-JpOmegaphi.x,-xsiOmega^2(Jt-Jp)cos(Omegat)\n Fx=mx..c-mepsilonOmega^2cos(Omegat+alfa)\n Fy=my..c-mepsilonOmega^2sin(Omegat+alfa)\n **RootLociJefDamp**\n Im(s)=wRe\n Re(s)=wIm\n Forward can be unstable\n Gup wFWDup wBWDdown\n Centrifugal stiffening wFWDeBWD up /(\n ceq={p+1 2}k,,Omege^(p-1)\n p=1viscous ')
    if (formulaSheet_chosen == 10):
        print(' **Reduction Tecniques**\n **Guyan**\n Mcond=M11+(K12*K22_1)*m22*(K22_1*K12^T)\n Kcond=K11-K12*K22_1*K21\n w=sqrt(Kcond/Mcond)\n **Errors2ModRed**\n -high freq modes neglected(if dof lot err=0)\n -modal coupling\n -if G&H!0 DONT\n **m!=n**\n multiply times phi*^T*M both sides\n Inertia forces related to slave degrees of freedom are actually not neglected, but their contribution to the kinetic energy is computed from a deformed configuration obtained on the basis of the master degrees of freedom alone.\n  ')

    if (formulaSheet_chosen == 11):
        print(' **Torsional Vibrations**\n leq=2*c+0.8*b+3/4*(D^4-d^4)/(D,^4-d,^4)*a+3/2*(D^4-d^4)/(b*s^3)*r\n Keqtor=GIpol/leq\n **KeqandJeq**\n K2=tau^2*k2\n same for J follow the transmission direction\n Jeqall = Jeqcrankshaft+mr1*r2^2+(mr2+mpiston)*r^2*f1(teta)+J0f2(teta)\n f1teta=(sinteta+alfasin(2teta)/(2cosbeta)-betacosteta/cosgamma)^2\n f2teta=alfa^2*(costeta/cosgamme)^2\n Mgas+Minertial=M=Jeqphiz..+0.5(w+phi.)^2dJeq/dteta\n Mgas=pAr*sqrt(f1(teta))\n teta=omega/2t (quatrotempi) \n Jb=Jrod=J0+mr1*a^2+mr2*b^2 \n 4tempi:wn=nOmega/2\n 2tempi:wn=nOmega\n angolo=k*4pi/2/n_cilindri  (2se4tempi e1se2tempi=stroke)\n J=Ip*lenght*rho\n Jeq=J+r^2*A*l*rho\n web = rectancgle\n mr2=mr*a/l\n mr1=mr*b/l\n alfa<0.3 beta=0\n Kgear=Kgear*tau^2\n Jgear=Jgear*tau^2  ')
        extra = 0
        print('load more?[1-y][0-n]')
        extra = int(input())
        if extra == 1:
            print(' F=a0+sum(aicos(2pii/T)+sum(bisin(2pii/T)\n a0=1/Tint(F dt)\n in 4 stroke \n odd-> gas p\n even-> gas p + inertiatorque\n tau=wout/win opposite to the one of MD\n **DampingT**\n Md=k,Ar^2\n Ccrank=k,Ar^2    ')
    if (formulaSheet_chosen == 12):
        print(' **FRFz**\n Dynamic Compliance:x/f m/N 0+in-in\n DynamicStiffness:f/x N/m 0-in+in\n Mobility:x./f -in+in-in\n Mc Impedance:f/x. Ns/m 1-in1\n Inertance:x../f m/s^2*N -in+in0\n Dynamic Mass:f/x.. s^2*N/m +in-in0\n n=dof=number of peaks\n antireson = n-i (11 ->1, 12->2,13->3) \n **accelerometerResponse**\n z=y0*r^2/sqrt((1-r^2)^2+(2*zita*r)^2)\n tgphi=-2*zita*r/(1-r^2)\n **1dof**\n H=1/k/sqrt((1-r^2)^2+(2*zita*r)^2)\n tanphi=-2zita*r/(1-r^2)\n Hij=xj0/Fk0\n j:output\n k:input  ')

    if (formulaSheet_chosen == 13):
        print(
            ' **MathStuff**\n d(cos)=-sin\n d(sin)=cos \n int(sin)=-cos\n int(cos)=sin\n cos pari\n sin dispari\n int(fg,)=fg-int(f,g)\n **TICTACMet**\n SDI\n +x^2e^x\n -2xe^x\n __.\n **2x2eigenproblem**  \n w^4+(-d11-d22)w^2+(d11*d22-d12*d21)\n **VectorProduct**\n +i -j +k repeat')
    print('need something else?[1-yes][0-no]')
    stop = int(input())
    if (stop == 0):
        print('see you later ;)')
    else:
        print('thanks Master GoodLuck ;)')
