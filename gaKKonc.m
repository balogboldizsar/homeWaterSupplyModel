function gaKKonc
nvars = 3;    % Number of variables
LB = [0,0,0];   % Lower bound
UB = [50,10,10];  % Upper bound
% nvars = 1;    % Number of variables
% LB = 0;   % Lower bound
% UB = 1000;  % Upper bound
options = optimoptions(@ga,'PlotFcn',{@gaplotbestf}, ...
    'Display','iter');
options.PopulationSize = 20;
[K,fval] = ga(@modell,nvars,[],[],[],[],LB,UB,[],options);
fprintf('\n K=[%5.3f,%5.3f,%5.3f];\n fval=%5.3f\n',K(1),K(2),K(3),fval);
end

function costFunc=modell(K)
global p0 ro g
global lambda L1 D1 A1 L2 D2 A2 dzeta Hv psz
global eMout p30
global nn n np jgPolE nMax nMin
global i

%termeszeti allandok
g=9.81; ro=1000; p0=1e5;

%rendszer jellemzoi
%csovek jellemzoi
L1=50; L2=0; %L2=100;
D1=32e-3; D2=D1/4; %D2=(25)*10^-3;
A1=D1^2*pi/4; A2=(rand()+0.05)*D2^2*pi/4;
lambda=0.02; dzeta=1; %dzeta=3;
 
%kut jellemzoi
Hv=20; psz=p0+ro*g*(L1-Hv);

%szivattyu jellemzoi
Qjg=(0:10:120)/60e3; %l/min->m^3/s
pjg=[80,78,75,72.5,70,67,63,59,52.5,46,41,35,29]*ro*g; %vom->Pa
jgPolE=polyfit(Qjg,pjg,2);
%plot(Qjg,polyval(jgPolE,Qjg),Qjg,pjg);
nn=3000; n=000; nMax=nn*1.2; nMin=0; %rpm

%szimulacio beallitasai
tmax=2; tt=0; dt=0.01;
i=1; j=2;

%szabalyozas parameterei
p3Ref=p0+2e5; %bar->Pa
Tm=0.01; TcPrev=Tm; %mintaveteli ido

%szabalyozas parameterei

p30=p0;
Qprev=0; Qnew=0;
eMout=zeros(floor(tmax/dt)+1,1);
eMout(1)=1-(p30-p0)/p3Ref;
pout=zeros(15*floor(tmax/dt)+1,1);
np=3000;

while tt<tmax
    if TcPrev>=Tm %mintaveteles szabalyozas
        n=control(p3Ref,p30,K(1),K(2),K(3));
        np=n;
        TcPrev=0;
    end
    
    [~,Qnew]=ode113(@ODE_Pump,[tt tt+dt],Qprev);
    Qnew(Qnew<0)=0;
    
    Qprev=Qnew(end);
    p30=p0+(dzeta+lambda*L2/D2)*ro/2/A2^2*Qprev^2+ro*L2/A2*(Qnew(end)-Qprev)/dt;
    
    pout(j:j-1+length(Qnew))=p30*ones(size(Qnew));
    
    tt=tt+dt;
    TcPrev=TcPrev+dt;
    i=i+1;
    j=j+length(Qnew);
end

eSout=1-(p30-p0)/p3Ref;
overShoot=-min(eSout);
if(overShoot<0)
    overShoot=0;
end
costFunc=sum(abs(eSout))/length(eSout)+0.5*overShoot;
Qstac=Qnew(end);

fprintf('K=[%02.2f,%02.2f,%02.2f] -> ',K(1),K(2),K(3));
fprintf('n=%5.0f rpm, Qst=%3.0f l/min, pst=%6.0f Pa\n',n,Qstac*60e3,p30);
end

function nOut=control(p3Ref,p30,Kp,Ki,Kd)
global eMout np nn i nMax nMin p0
e=1-(p30-p0)/p3Ref;
eMout(i+1)=e;
% if i>10
%     summ=sum(eMout(i-10:i));
% else
%     summ=0;
% end
if np<nMax
    summ=sum(eMout);
else
    summ=0;
end

if i>2
    diff=e-eMout(i);
else
    diff=0;
end

nc=Kp*e+Ki*summ/10+Kd*diff;
n=np+nc*nn;

if (n>nMax)
    n=nMax;
end
if(n<nMin)
    n=nMin;
end

nOut=n;
end

function dQdt=ODE_Pump(~,y)
global p0 ro
global lambda L1 D1 A1 L2 D2 A2 dzeta psz
Q=y(1);
p2=p0+dpsz(Q);
K1=p2-psz-p0;
K2=ro/2*(lambda*L1/D1/A1^2+(dzeta+lambda*L2/D2)/A2^2);
K3=ro*L1/A1+ro*L2/A2;
dQdt=(K1-K2*Q^2)/K3;
end

function out = dpsz(Q)
global nn n jgPolE
jgPol(1)=jgPolE(1);
jgPol(2)=jgPolE(2)*(n/nn);
jgPol(3)=jgPolE(3)*(n/nn)^2;
if Q>=0
    out=polyval(jgPol,Q);
    if (out<0)
        out=0;
    end
else
    out=jgPol(3)+(-1e5)*Q;
end
end