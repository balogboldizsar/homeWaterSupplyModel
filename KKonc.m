function KKonc
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
A1=D1^2*pi/4; A2=1*D2^2*pi/4;
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
tmax=1; tt=0; dt=0.01;
i=1; j=2;

%szabalyozas parameterei
p3Ref=p0+3e5; %bar->Pa
tM=0.01; TcPrev=tM; %mintaveteli ido
np=0;
%allando fojtas
 K=[10,0,0];
%K(1)=1.1479; K(2)=0.1; K(3)=0.1;
% K=[1.0951,0.4255,7.4180];
% K =[31.3045,9.9974,6.9625];
% K =[1.4160,0.2500,6.9625];
% K =[5.4172,2.3496,1.8633];
%K=[17.073,0.437,1.067]; %100m
%K=[18.909,9.222,4.271]; %100m
%K=[12.455,9.107,0.753]; %100m, nMax=4800
%K=[55.311,8.594,9.469]; %0m, dzeta=100, nMax=300*1.2;
%K=[0.010,0.002,1.646]; %0m, D2=D1/4;
%K=[49.369,9.779,6.993];
%K=[0.682,0.502,0.098];
%K=[0.630,0.237,0.165];
%K=[2.391,0.812,5.402];
%K=[59.71,90.02,45.81]
%K=[23.686,1.665,7.017]; %dolg
%K=[2.903,1.605,0.315]; %ea

%valtozo fojtas
%K =   7883    4.8712    6.4657
%K=[15.7099,7.0399,4.4946];
%K=[21.8822.32,6.8926,8.8026];
% K =[2.1993,0.2061,8.1116];
% K =[10.3849,2.9095,1.9403];
% K =[5.7268,2.1814,8.0093];
%K=[5.292,2.578,2.279];
%K=[0.336,0.070,0.109]; %dolg
%K=[0.841,0.000,0.617]; %ea
%K=[0.342,0.000,5.385];
%K=[3.964,2.849,3.857];

% K=[5.554,3.233,3.345]; %tulloves is
% K=[6.260,1.488,2.565];

%K=[36.439,4.147,8.023];
%K=[38.332,1.397,8.581];

%eredmenymatrixok
p30=p0;
tout=zeros(15*floor(tmax/dt)+1,1);
Qout=zeros(15*floor(tmax/dt)+1,1);
pout=ones(15*floor(tmax/dt)+1,1)*p0;
nout=zeros(15*floor(tmax/dt)+1,1);
eMout=zeros(floor(tmax/dt)+1,1);
eMout(1)=1-(p30-p0)/p3Ref;
tMout=zeros(floor(tmax/dt)+1,1);

while tt<tmax
    if TcPrev>=tM %mintaveteles szabalyozas
        n=control(p3Ref,p30,K(1),K(2),K(3));
        TcPrev=0;
        tMout(i)=tt;
    end
    [t,Qnew]=ode113(@ODE_Pump,[tt tt+dt],Qout(j-1));
    Qnew(Qnew<0)=0;
    Qout(j:j-1+length(t))=Qnew;

    p30=p0+(dzeta+lambda*L2/D2)*ro/2/A2^2*Qnew(end)^2+...
        ro*L2/A2*(Qnew(end)-Qout(j-1+length(t)))/dt;
    
    tout(j:j-1+length(t))=t;
    nout(j:j-1+length(t))=n*ones(size(t));
    pout(j:j-1+length(t))=p30*ones(size(t));
    
    p30=pout(j+length(t)-1);
    tt=tt+dt;
    TcPrev=TcPrev+dt;
    np=n;
    i=i+1;
    j=j+length(t);
end
Qout=Qout(1:j-1);
tout=tout(1:j-1);
pout=pout(1:j-1);
nout=nout(1:j-1);
eSout=1-(pout-p0)/p3Ref;
plotFunc(tout,pout,Qout,nout,eSout,tMout,eMout);
%plotWp(tout,pout,Qout)

Qstac=Qout(end);

fprintf('\n n=%5.3f rpm, Qstac=%g l/min\n',n,Qstac*60e3);
end

function nOut=control(p3Ref,p30,Kp,Ki,Kd)
global eMout np nn i nMax nMin p0
e=1-(p30)/p3Ref;
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

function plotWp(tout,pout,Qout)
figure()
subplot(2,1,1)
plot(tout,pout/1e5); grid on;
xlabel('t, (s)'); ylabel('p(end), (bar)'); 
ax=gca();
ax.FontSize=14;
subplot(2,1,2)
plot(tout,Qout*60e3); grid on;
xlabel('t, (s)'); ylabel('Q(end), (l/perc)')
ax=gca();
ax.FontSize=14;
end

function plotFunc(tout,pout,Qout,nout,eSout,tMout,eMout)
% figure()
% subplot(2,2,1)
% plot(tout,pout/1e5); grid on;
% xlabel('t, (s)'); ylabel('p3, (bar)')
% ax=gca(); ax.FontSize=12;
% subplot(2,2,2)
% plot(tout,Qout*60e3); grid on;
% xlabel('t, (s)'); ylabel('Q, (l/perc)')
% ax=gca(); ax.FontSize=12;
% subplot(2,2,3)
% plot(tout,eSout); grid on;
% xlabel('t, (s)'); ylabel('e, (Pa/Pa)')
% ax=gca(); ax.FontSize=12;
% subplot(2,2,4)
% plot(tout,nout); grid on;
% xlabel('t, (s)'); ylabel('n, (f/perc)')
% ax=gca(); ax.FontSize=12;

%figure()
% subplot(1,2,1)
% plot(tout,pout/1e5, 'LineWidth',2); grid on;
% xlabel('t, (s)'); ylabel('p3, (bar)'); 
% ax=gca();
% ax.FontSize=14;
% subplot(1,2,2)
% plot(tout,Qout*60e3, 'LineWidth',2); grid on;
% xlabel('t, (s)'); ylabel('Q, (l/perc)')
% ax=gca();
% ax.FontSize=14;
% subplot(2,2,3)
% plot(tout,nout, 'LineWidth',2); grid on;
% xlabel('t, (s)'); ylabel('n, (f/perc)')
% ax=gca();
% ax.FontSize=14;

% figure();
% subplot(2,1,1)
% plot(tout,eSout), grid on
% xlim([0,2])
% xlabel('t (s)'), ylabel('eS (Pa/Pa)')
% subplot(2,1,2)
% plot(tMout,eMout), grid on
% xlabel('t (s)'), ylabel('eM (Pa/Pa)')

figure()
subplot(2,2,1)
plot(tout,pout/1e5, 'k', 'LineWidth',1.5); grid on;
xlabel('t, (s)'); ylabel('p3, (bar)'); 
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
plot(tout,nout,'k' , 'LineWidth',1.5); grid on;
xlabel('t, (s)'); ylabel('n, (f/perc)')
ax=gca();
ax.FontSize=20;
ax.XLim(1)=0;
end