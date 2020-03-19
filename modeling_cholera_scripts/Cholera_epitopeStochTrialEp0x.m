
function []=Cholera_epitopeStochTrial
 
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
  
fitnessH=[.95,0.9,.925,.975,1.1,1.15,1,1,1];
fitnessE=[.76,.65,.72,.81,1.5,1.75,1,1,1];
   
% fitnessH=[.95,0.65,.75,.85,1,1,1,1,1];
% fitnessE=[.03,.01,.02,.05,1,1,1,1,1];
% fitnessH=[.95,0.65,.75,.85,1.15,1.25,1,1,1];  
% fitnessE=[.03,.01,.02,.25,1.75,2.5,1,1,1];
% fitnessH=[.75,0.75,.75,.75,1.15,1.25,1,1,1];
% fitnessE=[.02,.02,.02,.02,1.75,2.5,1,1,1];
Indp=2.^(0:1:d+l+h-1);
 ts=1.44;
amp=0.05;
amp2=0.5;
amp3=0.5;
mu=1/(365*55);
muP=mu*ts
gamma=1/7;
%muv(end)=10*mu;
bi=0.25/2.35;
biP=bi/N
bw=0.5;
bwP=bw/N
alpha=1/365;

xi=1/100;
xir=1.525;
xir=1;
K=1;
N=3e6;
sigma=0.36;
zeta=1;
ts=1.44;
 KP=sigma/zeta*N
 xiP=zeta/sigma*xi
eps=(5e-5)*ts;
sigs=.98;
r=0;
r=2;
rP=r*xi
nu=3.5;
nuP=nu*xi
R0h=(bi)/(gamma+mu)
R0c=(bw*xir)/(nu*(gamma+mu))
R0g=r/nu
R0=1/2*(R0h+R0g+sqrt((R0h-R0g)^2+4*R0c))
extinctlim=2;
  
X=zeros((m+2)*2,1);
X(1)=1;
%X=[b/a, 26, 0, 0, 0,0,0 0, 0,.1,.1,0.1]';
X(2)=.00005;
%X(2^(d)+2)=50;
  
%X(m+2:m+d+1)=0.1*ones(d,1);
  
  
  
t_r=1;
%pop_c=0.0001;
  
L2=length(X);
  
  
%numreac1=2^numEpi;
  
%[A,B] = hypercube_imm(d)
  
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
%v1=A+eye(numreac1);
filename = '/Users/cxb0559/Dropbox/ImmuneAgeModel/CholeraData/haiti_rainfall_water_chl_cases_data.xlsx';
xlRange = 'A21:E145';
dataArray_v = xlsread(filename,xlRange);
I=1:2:145-20;
Tr=dataArray_v(I,1);
 Cases=dataArray_v(I,5);
 Precip=dataArray_v(I,2);
Tr=Tr+365*(1903+11/12-2010.75);
 
  
 a1 =       116.9 ;
       b1 =      0.6061 ;
       c1 =      0.1992 ;
       a2 =       47.68  ;
       b2 =       6.072  ;
       c2 =       2.198  ;
       a3 =       42.07 ;
       b3 =        1.16  ;
       c3 =       1.719  ;
       a4 =       19.89  ;
       b4 =       7.097 ;
       c4 =       1.161 ;
        
  
 per=3.75;
 cavg=75;
 cmax=164;
  
  
  
%epsilon=.03;
  
  
 
  
tau=2;
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
  
Am=zeros(size(Af));
  
 Xv=X(2:m+1);
 Xv(sum(Bf(:,[1:d+h+l-3,d+h+l]),2)==0)=.000005*ones(2^2,1);
 Xv([0,0,0,0,0,1,0,0,0]*Indp'+1)=.00005;
  Xv([0,0,0,0,1,0,0,0,0]*Indp'+1)=.00005;
    %Xv([0,0,0,1,0,0,0,0,1]*Indp'+1)=.0005;
  
 %Xv(sum(Bf(:,1:d+h+l-1),2)==0)=.005*ones(2^1,1);
 X(2:m+1)=Xv;
%  size(Xv)
%  size(X(m+3:end))
%  size(aE)
%  size(aH.*((bw)*X(m+3:2*m+2)+(bi)*X(2:m+1)))
%  size(xi.*(X(2:m+1)+X(m+3:end).*(aE.*(ones(m,1)-X(m+3:end)/K)-ones(m,1))))
  
%%
  
options=odeset('AbsTol', 1e-25*ones(1,2*m+4));
  
function [dxdt] = myodea(t,x,bw,bi)
%    xt=mod(t/365,per);
%     c =  ((a1*sin(b1*xt+c1) + a2*sin(b2*xt+c2) + a3*sin(b3*xt+c3) + a4*sin(b4*xt+c4))-cavg)/cmax;
   c=interp1(Tr,Precip,t)-mean(Precip);
    c=c./max(abs(Precip-mean(Precip)));
dxdt(1) = (mu - mu*x(1)+3*alpha*x(m+2) - aH'*((bw*(1+amp*c))*x(m+3:2*m+2)+(bi)*x(2:m+1))*x(1))*ts;
dxdt(2:m+1) = (x(2:m+1).*(aH*bi*x(1) - gamma - mu) +(aH.*bw*(1+amp*c))*x(1).*x(m+3:2*m+2))*ts;
dxdt(m+2) = ((3*alpha)*x(2*m+4)-(mu+3*alpha)*x(m+2))*ts;
dxdt(m+3:2*m+2) = (xi.*(xir*x(2:m+1)*(1+amp2*c)+x(m+3:2*m+2).*(r*ones(m,1)*(1-sum(x(m+3:2*m+2))/K)-nu*(1-amp3*c)*aE.*ones(m,1))))*ts;
dxdt(2*m+3)=(sigs*gamma'*ones(1,m)*x(2:m+1)-(mu+3*alpha)*x(2*m+3))*ts;
dxdt(2*m+4)=((3*alpha)*x(2*m+3)-(mu+3*alpha)*x(2*m+4))*ts;
dxdt=dxdt';
end
  
  
%% 
%options=odeset('AbsTol', [1e-23,1e-23,1e-23,1e-25,1e-25,1e-25,1e-25]);
t=0;
tend=1882;
T=0;
Y=X;
Mutsum=zeros(m,1);
neH=zeros(d+l+h,1);
NeH=neH;
pepH=zeros(d+l+h,1);
PepH=pepH;
ne=zeros(d+l+h,1);
Ne=ne;
pep=zeros(d+l+h,1);
Pep=pep;
vlH=sum(X(2:m+1));
vl=vlH+sum(X(m+3:end));
while t<tend
%A1=eps*a1.*[X(1)*X(2);X(1)*X(3);X(1)*X(4);X(1)*X(5)];
  
                  
  
                    Tspan=t:1:t+tau;
                    if t==600
                        bw=.485*bw;
                        bi=.485*bi;
                    elseif t==900
                            bw=.58*bw;
                            bi=.58*bi;
                    end
                      
                    [T1, Z] = ode45(@(t,x)myodea(t,x,bw,bi),Tspan,X,options);
                    X=(Z(2:length(T1),:))';
                    T2=T1(2:end)';
                    Xv=N*X(2:m+1,:);
                    Ie=find(Xv<extinctlim);
                    Xv(Ie)=zeros(size(Xv(Ie)));
                    %Xv;
                    Xe=sigma*N/zeta*X(m+3:2*m+2,:);
                    Ie=find(Xe<extinctlim);
                    Xe(Ie)=zeros(size(Xe(Ie)));
%                      cv =   ((a1*sin(b1*mod(T2,1380)/7+c1) + a2*sin(b2*mod(T2,1380)/7+c2))-40)/70;
                    k=floor(sum(Xv,2));
                      
                    Xv=Xv(:,end);
                    X=X(:,end);
                    %k=poissrnd(tau*A1);
                    
                    %epsV=binornd(repmat(k,1,d),eps*ones(m,d));
                    
                      
                      
                    for i=1:m
                        epsV=binornd(k(i),ones(1,d+h+l)*eps);
                          
                        Im=find(Af(i,:)==1);
                        Am(i,Im)=Af(i,Im).*epsV;
                          
                          
                        Mutsum(i)=sum(epsV);
                    end
                         
                      
                      
                    Xv=max([Xv-Mutsum,zeros(size(Xv))],[],2);
                    Xv=Xv+(ones(1,m)*Am)';
                    VLH=sum(Xv);
%                     if VLH==0
%                         break;
%                     end
                    for j=1:d+l+h
                       p=sum(Xv(Bf(:,j)==0))/VLH;
                        PepH(j)=p;
                         
                 %G=1-2*p*(1-p);
                    pi=2*p*(1-p);
                    %theta=1/G-1;
                    %Ne(j)=theta/(2*(2*eps));
                    NeH(j)=pi/(2*eps);
                          
                    end
                      
%                    Xe=Xe(:,end);
                    %Xv;
%                     size(ones(m,1)*(1-zeta/(N*sigma)*sum(Xe(1:m))/K))
%                     size(xi*aE.*Xe(1:m))
                      %k=floor(xi*Xe(1:m).*ones(m,1)*(1-zeta/(N*sigma)*sum(Xe(1:m))/K));
                        
                    k=floor(sum(xi*r*Xe.*ones(m,length(T1)-1).*repmat(ones(1,length(T1)-1)-zeta/(N*sigma)*sum(Xe)/K,m,1),2));
                      
                    %k=poissrnd(tau*A1);
                      
                    %epsV=binornd(repmat(k,1,d),eps*ones(m,d));
                    
                      Xe=Xe(:,end);
                      
                    for i=1:m
                        epsV=binornd(k(i),ones(1,d+h+l)*eps);
                          
                        Im=find(Af(i,:)==1);
                        Am(i,Im)=Af(i,Im).*epsV;
                          
                          
                        Mutsum(i)=sum(epsV);
                    end
                         
                      
                      
                    Xe=max([Xe-Mutsum,zeros(size(Xe))],[],2);
                    Xe=Xe+(ones(1,m)*Am)';
                    VL=sum(Xe)+VLH;
                    if VL==0
                        break;
                    end
                    for j=1:d+l+h
                       p=sum(Xv(Bf(:,j)==0)+Xe(Bf(:,j)==0))/VL;
                        Pep(j)=p;
                         
                 %G=1-2*p*(1-p);
                    pi=2*p*(1-p);
                    %theta=1/G-1;
                    %Ne(j)=theta/(2*(2*eps));
                    Ne(j)=pi/(2*eps);
                          
                    end
                      
                      
                    X(2:m+1)=Xv/N;
                    X(m+3:2*m+2)=zeta/(N*sigma)*Xe;
                      
                      
                    t=t+tau;
                      
               % A1
                %sum(A1)
               %r(2)*sum(A1)
                  
                     
            
            
          
          
Y=[Y,X];
T=[T,t]; 
vl=[vl,VL];
ne=[ne,Ne];
pep=[pep,Pep];
vlH=[vlH,VLH];
neH=[neH,NeH];
pepH=[pepH,PepH];
  
if isnan(X)~=zeros(size(X))
    break
end
end  
 save('timeD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat','T')
save('stYD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat','Y')
save('neD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat','ne')
save('pepD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat','pep')
save('neHD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat','neH')
save('pepHD1tdccsn2Ahcs2fAsEc2copy2Epcs23fc2cb.mat','pepH')
figure(1)
plot(T,vlH/4)
  
for i=2:2^(d+h+l)+1
figure(2)
plot(T,N*Y(i,:),'LineWidth',3);
ylim([0 20000]);
hold on
end
legend(BinS(1:2^(d+h+l)))
%legend('show')
figure(3)
plot(T,vl)
  
  
figure(4)
plot(T,mean(ne),'LineWidth',3)
  
figure(5)
plot(T,mean(ne(d+h+1:end,:)),'LineWidth',3)
%figure(5)
%plot(T,ne(d+l,:))
for i=1:l+d+h
figure(6)
plot(T,pep(i,:),'LineWidth',3)
%ylim([0 20000])
hold on
  
end
  
for i=1:l+d+h
figure(7)
plot(T,pepH(i,:),'LineWidth',3)
%ylim([0 20000])
hold on
  
end
for i=m+3:2*m+2
figure(8)
plot(T,N*Y(i,:),'LineWidth',3);
ylim([0 20000]);
hold on
end
legend(BinS(1:2^(d+h+l)))
 
Yh=Y(2:m+1,:);
 
Ye=Y(m+3:2*m+2,:);
mYh=mean(Yh,2);
 
mYe=mean(Ye,2);
 
[Ms,Is]=sort(Yh,'descend');
% mYh(1:20)
[Mmh,Imh]=sort(mYh,'descend');
 
% Imh(1:20)
% m1=Mmh(1:20);
Ist=unique(Is(1,:));
% m2=max(Yh(Imh(1:20),:),2);
for i=1:length(Ist)
    figure(11)
    plot(T,N*Yh(Ist(i),:),'LineWidth',3);
    hold on 
end
legend(BinS(Ist));
 
[Mme,Ime]=sort(mYe,'descend');
 
for i=1:20
figure(9)
plot(T,N*Yh(Imh(i),:),'LineWidth',3);
ylim([0 1000]);
hold on
end
legend(BinS(Imh(1:20)))
 
for i=1:20
figure(10)
plot(T,sigma*N/zeta*Ye(Ime(i),:),'LineWidth',3);
ylim([0 1000]);
hold on
end
legend(BinS(Ime(1:20)))
end
      
        
    
       
      
      
     
      
        
    
       
      
     