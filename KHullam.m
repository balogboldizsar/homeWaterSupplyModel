function KHullam
global p0 ro g
global lambda L1 D1 A1 Af Hv psz
global i elemSzam B a dt dx pPrev vPrev Q eMout pRef
global nn n np jgPolE nMax nMin
global options

options = optimoptions('fsolve','Display','none');

%termeszeti allandok
g=9.81; ro=1000; p0=1e5;

%rendszer jellemzoi
%csovek jellemzoi
L1=50; L2=0; %L2=100; 
D1=32e-3; D2=D1/4; %D2=(25)*10^-3; 
A1=D1^2*pi/4; Af=0.25*D2^2*pi/4; %Af=D2^2*pi/4; 
lambda=0.02; dzeta=3; %dzeta=1;

%kut jellemzoi
Hv=20; psz=p0+ro*g*(L1-Hv);

%szivattyu jellemzoi
Qjg=(0:10:120)/60e3; %l/min->m^3/s
pjg=[80,78,75,72.5,70,67,63,59,52.5,46,41,35,29]*ro*g; %vom->Pa
jgPolE=polyfit(Qjg,pjg,2);
%plot(Qjg,polyval(jgPolE,Qjg),Qjg,pjg);
nn=3000; n=000; nMax=nn*1.2; nMin=0; %rpm

%hullamterjedes parameterei
B=2.1e9;
a=sqrt(B/ro);
elemSzam=20;
dx=L1/elemSzam;
dt=dx/a;

%szimulacio parameterei
tmax=2; t=-0.5;
i=1;

%szabalyozas parameterei
pRef=p0+3e5;
tM=0.01; TcPrev=tM; %mintaveteli ido
np=0;

%K=[0.01,0,0];
%K=[0.072,0.001,0.713]; %fval=0.117
%K=[1.617,8.697,7.208]; %valtozo fojtas
%K=[0.101,9.339,2.162]; %valtozo fojtas
K=[0.095,1.022,4.036]; %opt4


%eredmenymatrixok
pPrev=((ro*g*L1):-(ro*g*L1)/elemSzam:0)+p0;
vPrev=zeros(elemSzam,1);
Q=zeros(elemSzam,1);
vnew=zeros(1,elemSzam);
pnew=vnew;
eMout=zeros(1,floor(tmax/dt)+1);
nout=zeros(1,floor(tmax/dt));

tout=zeros(floor(tmax/dt)+1,1);
pout=zeros(elemSzam,floor(tmax/dt)+1);
Qout=zeros(elemSzam,floor(tmax/dt)+1);

while t<tmax
    if TcPrev>=tM&&t>0 %mintaveteles szabalyozas
        n=control(pRef,pPrev(end),K(1),K(2),K(3));
        np=n;
        TcPrev=0;
    end
    
    [pnew(1),vnew(1)]=BCLeft();
    [pnew(2:elemSzam-1),vnew(2:elemSzam-1)]= Update();
    [pnew(end),vnew(end)]=BCRight();
    
    pPrev=pnew;
    vPrev=vnew;
    
    pout(:,i)=pnew;
    Qout(:,i)=A1*vnew;
    tout(i)=t;
    nout(i)=n;
    t=t+dt;
    TcPrev=TcPrev+dt;
    i=i+1;
end
eSout=1-pout(end,:)/pRef;

plotFunc(tout,pout,Qout,eSout,nout);
end


function nOut=control(pRef,pPrev,Kp,Ki,Kd)
global eMout i np nn nMax nMin
e=1-pPrev(end)/pRef;
eMout(i)=e;
% if i>10
%     summ=sum(eMout(i-10:i));
% else
%     summ=0;
% end
if np<nMax
    summ=sum(eMout)/length(eMout);
else
    summ=0;
end

if i>2
    diff=e-eMout(i);
else
    diff=0;
end

nc=Kp*e+Ki*summ+Kd*diff;
n=np+nc*nn;

if (n>nMax)
    n=nMax;
end
if(n<nMin)
    n=nMin;
end

nOut=n;
end

function [pnewL,vnewL]=BCLeft()
global ro a psz A1 Kl pPrev vPrev
global do_runtime_plot options

Kl=pPrev(2)-ro*a*vPrev(2);
Qsz=fsolve(@eqBCLeft,A1*vPrev(1),options);

vnewL=Qsz/A1;
if vnewL(1)<0 %visszacsapo szelep
    vnewL(1)=0;
end
pnewL=Kl+ro*a*vnewL;

if do_runtime_plot==1
    plot_jg(pnewL-psz,Qsz);
end
end

function out = eqBCLeft(Q)
global Kl ro lambda a D1 A1 psz g dt vPrev
S_betaR=-(-g-lambda/2/D1*abs(vPrev(2))*vPrev(2));
out=Kl+ro*a*Q/A1-(psz+dpsz(Q))+dt*a*ro*S_betaR;
end

function [pnew,vnew]=Update()
global elemSzam D1 dt ro a g lambda pPrev vPrev
pnew=zeros(elemSzam-2,1);
vnew=pnew;
for j=2:(elemSzam-1)
    S_alpha=-g-lambda/2/D1*abs(vPrev(j-1))*vPrev(j-1);
    S_beta =-(-g-lambda/2/D1*abs(vPrev(j+1))*vPrev(j+1));
    alpha=pPrev(j-1)+ro*a*vPrev(j-1)+dt*ro*a*S_alpha;
    beta =pPrev(j+1)-ro*a*vPrev(j+1)+dt*ro*a*S_beta;
    pnew(j-1)=(alpha+beta)/2;
    vnew(j-1)=(alpha-beta)/2/ro/a;
end
end

function [pnewR,vnewR]=BCRight()
global Af ro p0 A1 pPrev options g vPrev lambda D1 a dt

% %eloirt nyomas
% S_alphaL=-g-lambda/2/D1*abs(vPrev(end-1))*vPrev(end-1);
% pnewR=p0+0*2e5;
% vnewR=(pPrev(end-1)+ro*a*vPrev(end-1)+dt*ro*a*S_alphaL-pnewR)/ro/a;

%eloirt aramlasi sebesseg
% S_alphaL=-g-lambda/2/D1*abs(vPrev(end-1))*vPrev(end-1);
% vnewR=50/(60e3)/A1;
% pnewR=pPrev(end-1)+ro*a*vPrev(end-1)+dt*S_alphaL-ro*a*vnewR;

%eloirt fojtas
pnewR=fsolve(@eqBCRight,pPrev(end),options);
dp=pnewR-p0;
vnewR=(Af*sqrt(2/ro*abs(dp)))*sign(dp)/A1;
end

function out = eqBCRight(pnewR)
global dt ro g a D1 lambda p0 A1 Af pPrev vPrev
dp=pnewR-p0;
vnewR=(Af*sqrt(2/ro*abs(dp)))*sign(dp)/A1;
S=dt*ro*a*(-g-lambda/2/D1*abs(vnewR)*vnewR);    %source term
out=(pnewR+ro*a*vnewR)-(pPrev(end-1)+ro*a*vPrev(end-1))-S;
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

function plotFunc(tout,pout,Qout,eout,nout)
figure()
subplot(2,2,1)
plot(tout,pout(end,:)/1e5, 'k', 'LineWidth',1.5); grid on;
xlabel('t, (s)'); ylabel('p, (bar)'); 
ax=gca();
ax.FontSize=20;
ax.XLim(1)=0;
% subplot(2,1,2)
% plot(tout,Qout(end,:)*60e3, 'LineWidth',2); grid on;
% xlabel('t, (s)'); ylabel('Q(end), (l/perc)')
% ax=gca();
% ax.FontSize=14;
% ax.XLim(1)=0;
subplot(2,2,2)
plot(tout,nout(end,:),'k' , 'LineWidth',1.5); grid on;
xlabel('t, (s)'); ylabel('n, (f/perc)')
ax=gca();
ax.FontSize=20;
ax.XLim(1)=0;

% figure();
% subplot(2,2,1)
% plot(tout,pout(end,:)/1e5); grid on;
% xlabel('t, (s)'); ylabel('p(end), (bar)');
% subplot(2,2,2)
% plot(tout,Qout(end,:)*60e3); grid on;
% xlabel('t, (s)'); ylabel('Q(end), (l/min)');
% subplot(2,2,3)
% plot(tout,eout); grid on;
% xlabel('t, (s)'); ylabel('e, (%)');
% subplot(2,2,4)
% plot(tout,nout); grid on;
% xlabel('t, (s)'); ylabel('n, (rpm)');
end