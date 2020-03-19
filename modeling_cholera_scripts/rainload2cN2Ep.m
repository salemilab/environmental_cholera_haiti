numEnv=4;
numhost=2;
neutAlleles=3;
q=4; % numEnv & neutral for host env.
d=numEnv;
l=neutAlleles;
h=numhost;
m=2^(d+h+l);
%fitnesscost=unifrnd(0.9,1,1,numEpi);
%fitnesscost=0.975*ones(1,numEpi);
% fitnessH=[1,1,1,1,0.85,1.01,1.5,1];
% fitnessE=[1.15,1.15,1.15,1.15,1.15,.85,1,1,1];
 
% fitnessH=[.95,0.9,.92,.97,1.1,1.15,1,1,1];
% fitnessE=[.75,.65,.67,.8,1.5,1.75,1,1,1];
  
fitnessH=[.95,0.9,.925,.975,1.1,1.15,1,1,1];
fitnessE=[.75,.65,.7,.8,1.5,1.75,1,1,1];

N=3e6;
sigma=0.36;
amp=0.05;
amp2=0.5;
amp3=0.05;
zeta=1;
ts=1/2;
filename = '/Users/cxb0559/Dropbox/ImmuneAgeModel/CholeraData/skyride_plot_data.xlsx';
xlRange = 'A3:E102';
dataArray = xlsread(filename,xlRange);

Te=flipud(dataArray(:,1));
Ne=flipud(dataArray(:,2));
NeL=flipud(dataArray(:,5));
NeU=flipud(dataArray(:,4));
eps=(5e-5)*ts;
filename = '/Users/cxb0559/Dropbox/ImmuneAgeModel/CholeraData/haiti_rainfall_water_chl_cases_data.xlsx';
xlRange = 'A21:E145';
dataArray_v = xlsread(filename,xlRange);
I=1:2:145-20;
Tdat=dataArray_v(I,1);
 Tr=Tdat+365*(1903+11/12-2010.75);
Tdat=Tdat/365+1903+11/12;
 Precip=dataArray_v(I,2);
 Casesd=dataArray_v(I,5)
 
[Af,Bf] = hypercube_imm(d+h+l);
 
A=Af(:,1:d+h);
B=Bf(:,1:d+h);
Bc=1-B;
Bcf=1-Bf;
Bin=num2str(Bf);
BinS= cell(m,1);
for i=1:m
    BinS{i,1}=strcat(Bin(i,:));
end
a1e=ones(size(Bf));
a1h=ones(size(Bf));
%tt=[fitnesscostM(:,k).*B(:,1),Bc(:,1)]
for k=1:d+h
a1e(:,k)=max([fitnessE(k)*B(:,k),Bc(:,k)],[],2);
a1h(:,k)=max([fitnessH(k)*B(:,k),Bc(:,k)],[],2);
%ind0=find(Bf(:,k)==0);
%X(ind0)=(b-a*X(1))/(beta0*X(1))/length(ind0);
 
 
%X(ind0)=5;
end
aH=(prod(a1h,2)).^(.65);
aE=(prod(a1e,2)).^(.65);
Td=load('timeD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat');
Yd=load('stYD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat');
ne1=load('neD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat');
pep1=load('pepD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat');

Td0=load('timeD1tdccsn2Ahcs2fAsEc2copy2Epcs23fcb0xc3c.mat');
Yd0=load('stYD1tdccsn2Ahcs2fAsEc2copy2Epcs23fcb0xc3c.mat');
ne0=load('neD1tdccsn2Ahcs2fAsEc2copy2Epcs23fcb0xc3c.mat');

% Td0n=load('timeD1tdccsn2Ahcs2fAsEc2copy20xn.mat');
% Yd0n=load('stYD1tdccsn2Ahcs2fAsEc2copy20xn.mat');
% ne0n=load('neD1tdccsn2Ahcs2fAsEc2copy20xn.mat');
%pep0=load('pepD1tdccsn2A0.mat');
T1=Td.T;
T00=Td0.T;
Y1=Yd.Y;

ne=ne1.ne/1000/(3/4.5)*2;

 c=interp1(Tr,Precip,T1)-mean(Precip);
    c=c./max(abs(Precip-mean(Precip)));

nez=ne0.ne/1000/(3/4.5)*2;
nez0=ne0n.ne/1000/(3/3.5)*2;
pep=pep1.pep;
Yh=N*Y1(2:m+1,:);
Yss=Y1(1,:);
Ys=N*Yss;
Cases=sum(Yh);
Yr=N*(Y1(m+2,:)+Y1(2*m+3,:)+Y1(2*m+4,:));

mu=1/(365*55);
gamma=1/7;
%muv(end)=10*mu;
bi=0.25/2.35;
bw=0.5;
biv=bi*ones(size(Yh));
bwv=bw*ones(size(Yh));
biv(:,T1>=600)=.485*biv(:,T1>=600);
bwv(:,T1>=600)=.485*bwv(:,T1>=600);
biv(:,T1>=900)=.58*biv(:,T1>=900);
bwv(:,T1>=900)=.58*bwv(:,T1>=900);
alpha=1/365;
xi=1/100;
K=1;
N=3e6;
sigma=0.36;
zeta=1;
ts=1;
Yw=sigma*N/zeta*Y1(m+3:2*m+2,:);
Yws=sum(Yw);
eps=(5e-5)*ts;
sigs=.98;
r=2;
%r=0;
nu=3.5;
R0h=(bi)/(gamma+mu)
R0c=(bw)/(nu*(gamma+mu))
R0g=r/nu
Km=[R0h,bw/nu;1/(nu*(gamma+mu)),R0g]

R0=1/2*(R0h+R0g+sqrt((R0h-R0g)^2+4*R0c))
% r11=(Km(1,1)+Km(1,2)*Km(2,1)*nu/(gamma+mu+nu))*1/gamma;
% r12=(Km(1,2)+Km(1,2)*Km(2,1)*(gamma+mu)/(gamma+mu+nu)+Km(1,2)*Km(2,2)*2)*1/nu;
size(aH')
size(bi*Yh)
size(Yss)
r11=aH'*(biv.*Yh).*Yss+xi*Cases.*(1+amp2*c);
r12=(1+amp*c).*(aH'*(bwv.*Yw)).*Yss+xi*r*Yws.*(1-zeta/N/sigma*Yws./K);
p11=r11./(r11+r12)
p12=r12./(r11+r12)
lambda=2*((p11.^2.*(aH'*(biv.*Yh))./Cases.^2+p11.*p12.*((1+amp*c).*(aH'*(bwv.*Yw)))./(Cases.*Yws)).*Yss+(p12.^2*xi*r)./Yws.*(1-zeta/N/sigma*Yws./K));
NeNa=1./lambda*(1e-5);
figure(15)
plot(T1,NeNa)

Cases=Cases/4;
Yh(1:20,:);
% % figure(10)
% % plot(T1,Cases./(.01*Ys))
% NeNa=Cases./(.01*Ys);
fn=.15;

%Ye=Y1(m+3:2*m+2,:);
mYh=mean(Yh,2);
[mah,Imh]=max(Yh);
It=unique(Imh,'stable');
T1=T1/365+2010.75;
T0=T00/365+2010.75;
%mYe=mean(Ye,2);
% mYh(1:20)
[Mmh,Imh]=sort(mYh,'descend');

% Imh(1:20)
% m1=Mmh(1:20);
% 
% m2=max(Yh(Imh(1:20),:),2);

%[Mme,Ime]=sort(mYe,'descend');

nem=mean(ne);
nem0=mean(nez);
NeT=nem*(1-fn)+NeNa*fn;
Tabm=[T1',nem',Cases'];
Tab=array2table(Tabm,...
    'VariableNames',{'Time','Ne','Cases'});
figure(11)
plot(T1,nem,'LineWidth',3)

figure(12)
plot(T1,Ys,T1,Yr)
figure(13)
plot(T1,NeT)

figure(1)
plotyy(T1,Cases,T1,nem)

figure(2)
plot(T1,pep)
legend;

figure(3)
plot(Te,Ne,T1,nem,T0,nem0)
hold on
jbfill(Te,NeU,NeL)

% figure(3)
% plot(Te,Ne,T1,nem)
% hold on
% jbfill(Te,NeU,NeL)

figure(4)
plot(Tdat,Casesd,T1,Cases)

for i=1:20
figure(9)
plot(T1,Yh(Imh(i),:),'LineWidth',3);
%ylim([0 1000]);
hold on
end
legend(BinS(Imh(1:20)))

figure(10)
plot(T1,Yh(It,:),'LineWidth',3);
%ylim([0 1000]);
legend(BinS(It))

cases7(Tdat,Casesd,T1,Cases)


NeCasesfig(T1,Cases,nem)
Y1=Yd0.Y;



 c=interp1(Tr,Precip,T00)-mean(Precip);
    c=c./max(abs(Precip-mean(Precip)));


Yh=N*Y1(2:m+1,:);
Yss=Y1(1,:);
Ys=N*Yss;
Cases=sum(Yh);
biv=bi*ones(size(Yh));
bwv=bw*ones(size(Yh));
biv(:,T0>=600)=.55*biv(:,T0>=600);
bwv(:,T0>=600)=.55*bwv(:,T0>=600);
% biv(:,T0>=900)=.8*biv(:,T0>=900);
% bwv(:,T0>=900)=.8*bwv(:,T0>=900);

mu=1/(365*55);
gamma=1/7;
%muv(end)=10*mu;
% bi=0.25/2.25;
% bw=0.5;
% alpha=1/365;
% xi=1/100;
% K=1;
% N=3e6;
% sigma=0.36;
% zeta=1;
% ts=1;
Yw=sigma*N/zeta*Y1(m+3:2*m+2,:);
Yws=sum(Yw);
eps=(5e-5)*ts;
sigs=.98;
r=0;

nu=2;
R0h=(bi)/(gamma+mu)
R0c=(bw)/(nu*(gamma+mu))
R0g=r/nu
Km=[R0h,bw/nu;1/(nu*(gamma+mu)),R0g]

R0=1/2*(R0h+R0g+sqrt((R0h-R0g)^2+4*R0c))
% r11=(Km(1,1)+Km(1,2)*Km(2,1)*nu/(gamma+mu+nu))*1/gamma;
% r12=(Km(1,2)+Km(1,2)*Km(2,1)*(gamma+mu)/(gamma+mu+nu)+Km(1,2)*Km(2,2)*2)*1/nu;
size(aH')
size(bi*Yh)
size(Yss)
xir=1.5244;
r11=aH'*(biv.*Yh).*Yss+xi*xir*Cases.*(1+amp2*c);
r12=(1+amp*c).*(aH'*(bwv.*Yw)).*Yss+xi*r*Yws.*(1-zeta/N/sigma*Yws./K);
p11=r11./(r11+r12)
p12=r12./(r11+r12)
lambda=2*((p11.^2.*(aH'*(biv.*Yh))./Cases.^2+p11.*p12.*((1+amp*c).*(aH'*(bwv.*Yw)))./(Cases.*Yws)).*Yss+(p12.^2*xi*r)./Yws.*(1-zeta/N/sigma*Yws./K));
NeNa=1./lambda*(1e-5);
NeT0=nem0*(1-fn)+NeNa*fn;


Ne2fig(T1,NeT,Te,Ne)
TeU=flipud(Te);
TeP=[Te;TeU];
NeP=[NeU;flipud(NeL)];
Ne3figE(Te,Ne,T1,NeT,T0,NeT0,NeP,TeP)
% hold on
% jbfill(Te,NeU,NeL)
% hold off
% for i=1:20
% figure(10)
% plot(T1,sigma*N/zeta*Ye(Ime(i),:),'LineWidth',3);
% ylim([0 1000]);
% hold on
% end
% legend(BinS(Ime(1:20)))

