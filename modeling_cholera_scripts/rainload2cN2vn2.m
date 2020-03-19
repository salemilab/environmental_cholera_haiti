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
  ss=(1882-1606)/2;
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
aH=prod(a1h,2);
aE=prod(a1e,2);
Td=load('timeD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat');
Yd=load('stYD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat');
ne1=load('neD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat');
pep1=load('pepD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat');

Td0=load('timeD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv04.mat');
Yd0=load('stYD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv04.mat');
ne0=load('neD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv04.mat');

Td0c=load('timeD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv01e1.mat');
Yd0c=load('stYD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv01e1.mat');
ne0c=load('neD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv01e1.mat');

Td0n=load('timeD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv01.mat');
Yd0n=load('stYD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv01.mat');
ne0n=load('neD1tdccsn2Ahcs2fAsEc2copy2nEoepcbv01.mat');
%pep0=load('pepD1tdccsn2A0.mat');
T1=Td.T;
T00=Td0.T;
T2=Td0n.T;
Y1=Yd.Y;
Y2=Yd0n.Y;
Y3=Yd0.Y;
T3=Td0.T;
Y4=Yd0c.Y;
T4=Td0c.T;
ne=ne1.ne/1000/(3/3.5)*2;

 c=interp1(Tr,Precip,T1)-mean(Precip);
    c=c./max(abs(Precip-mean(Precip)));

nez=ne0.ne/1000/(3/3.5)*2;
nez0=ne0n.ne/1000/(3/3.5)*2;
pep=pep1.pep;
Yh=N*Y1(2:m+1,:);
%T1=T1(1:end-ss);
Yh2=N*Y2(2:m+1,2:end);
Y2s=N*Y2(1,:);
Yh3=N*Y3(2:m+1,2:end);
Y3s=N*Y3(1,:);
Y4s=Y4(1,:)*N;
Yh4=N*Y4(2:m+1,2:end);
T2=T2(2:end);
T3=T3(2:end);
T4=T4(2:end);
Yss=Y1(1,:);
Ys=N*Yss;
Cases=sum(Yh);
Cases2=sum(Yh2)/4;
Cases3=sum(Yh3)/4;
Cases4=sum(Yh4)/4;

mu=1/(365*55);
gamma=1/7;
%muv(end)=10*mu;
bi=0.25/2.25;
bw=0.5;
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
nu=3.5;
R0h=(bi)/(gamma+mu)
R0c=(bw)/(nu*(gamma+mu))
R0g=r/nu
Km=[R0h,bw/nu;1/(nu*(gamma+mu)),R0g]

R0=1/2*(R0h+R0g+sqrt((R0h-R0g)^2+4*R0c))
% r11=(Km(1,1)+Km(1,2)*Km(2,1)*nu/(gamma+mu+nu))*1/gamma;
% r12=(Km(1,2)+Km(1,2)*Km(2,1)*(gamma+mu)/(gamma+mu+nu)+Km(1,2)*Km(2,2)*2)*1/nu;


Cases=Cases/4;
Yh(1:20,:);
figure(10)
plot(T1,Cases./(.01*Ys))
NeNa=Cases./(.01*Ys);
fn=.15;

Ye=Y1(m+3:2*m+2,:);
mYh=mean(Yh,2);
T1=T1/365+2010.75;
T2=T2/365+2010.75;
T3=T3/365+2010.75;
T4=T4/365+2010.75;
T0=T00/365+2010.75;
mYe=mean(Ye,2);
mYh(1:20)
[Mmh,Imh]=sort(mYh,'descend');

Imh(1:20)
m1=Mmh(1:20);

m2=max(Yh(Imh(1:20),:),2);

[Mme,Ime]=sort(mYe,'descend');



figure (11)
plot(T0,Y3s)




% figure(3)
% plot(Te,Ne,T1,nem)
% hold on
% jbfill(Te,NeU,NeL)

% figure(4)
% plot(Tdat,Casesd,T1,Cases)



% cases7(Tdat,Casesd,T1,Cases)
% hold on
% plot(T2,Cases2,T3,Cases3,T4,Cases4)
% hold off
% 
% 
casesvac(Tdat,Casesd,T1,Cases,T2,[Cases3;Cases2],T4,Cases4)
for i=1:20
figure(9)
plot(T1,sigma*N/zeta*Ye(Ime(i),:),'LineWidth',3);
ylim([0 1000]);
hold on
end
legend(BinS(Ime(1:20)))

