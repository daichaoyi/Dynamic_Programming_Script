clear all;
close all;
tic
  
% global alpha beta lambda rho p y r ki ai T amat amin amax ks v t i j asta acon aconi ci c c1 s theta

% set parameters
alpha=0.5; 
beta=0.75;
T=10;              % life span
lambda=0.2;
rho=1;
p=1;               % price for one unit of consumption
y=3;               % income
r=0.05;            % interest rat
ki=5;             % initial k in period t=0
ai=0;              % initial a in period t=0

% Construct a Grid for initial state of Kss
kmin=23;            % minimum
kmax=24;           % maximum
kgrid=40;          % grid points + 1
grid1=(kmax-kmin)/kgrid;      % grid
kmat=kmin:grid1:kmax;
kmat=kmat';        % make column vector as opposed to row

% Construct a Grid for initial state of A 
amin=0;      
amax=50;
agrid=30;
grid2=(amax-amin)/agrid;      % grid
amat= amin:grid2:amax;
amat= amat';

[N, L]=size(kmat);
[M, U]=size(amat);

aconi=zeros(T+1,N);
ci=zeros(1,N);
vi=zeros(T+1,N);

parfor s=2:9
   asta = zeros(T+1,M);
   acon = zeros(T+1,N,M);
   c = zeros(T+1,N,M);
   v = zeros(T+1,N,M);            % set v0 as a n*m matrix
   theta = (1/(T-s))*1.05;
   t=T+1;
  while t>0
    for i=1:N
        ks=kmat(i,1);
        for j=1:M  
         if t>T                                                         % find optimal VT
            acon(t,i,j)= amat(j,1);                                     % try grid AT
            c(t,i,j)=((1+r)*acon(t,i,j)-theta*(1-lambda)*ks+y)/p;       % calculate optimal CT
            if c(t,i,j)<0
            v(t,i,j)=-999;
            else
            v(t,i,j)=(1/rho)*(c(t,i,j)^alpha+beta*(ki+ks)^alpha)^(rho/alpha);       % calculate the optimal 
            end
           
           elseif t>s+1                                 % find optimal As+2 to AT
             asta(t,j)= amat(j,1)                       % try grid A
             acon(t,i,j)=optsearch(1,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta);                        % find optimal At+1
             c(t,i,j)=((1+r)*asta(t,j)-acon(t,i,j)-theta*(1-lambda)*ks+y)/p;        % calculate optimal (Cs+1-CT-1)
             v(t,i,j)=utility1(alpha,beta,lambda,rho,p,y,r,ki,amat,T,acon(t,i,j),asta,ks,i,j,t,v,c,theta);            % calculate optimal aggregate utility from t    
 
           elseif t>s                                    % find optimal As+1
              asta(t,j)= amat(j,1)                       % try grid As 
              acon(t,i,j)=optsearch(2,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta);                              % find optimal As+1
              c(t,i,j)=((1+r)*asta(t,j)-acon(t,i,j)-lambda*ks+y)/p                % calculate optimal Cs
              v(t,i,j)=utility2(alpha,beta,lambda,rho,p,y,r,ki,amat,acon(t,i,j),asta,ks,i,j,t,v);  
  
              elseif t>1                           % find optimal A2 to As; V1
                 asta(t,j)= amat(j,1)                                               % try grid At
                 acon(t,i,j)= optsearch(3,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta)                          % find optimal At+1
                 c(t,i,j)= ((1+r)*asta(t,j)-acon(t,i,j)+y)/p;                       % calculate optimal Ct(C1-Cs-1)
                 v(t,i,j)= utility3(alpha,beta,rho,p,y,r,ki,amat,acon(t,i,j),asta,i,j,t,v);                                  % calculate optimal aggregate utility from t
 
              else
                 aconi(s,i)= optsearch(4,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta);        % find optimal A1; V0(vi)
                 ci(s,i)=((1+r)*ai-aconi(s,i)+y)/p;                                                                % find optimal C0   
                 vi(s,i)= utility4(alpha,beta,rho,p,y,r,ki,ai,amat,aconi(s,i),i,t,v);                  % the optimal value of aggregate utility form t=0
          end
        end
      end
    t=t-1;
  end
end


% find the optimal consumption and asset path
[num] = max(vi(:));
[opts,optks] = ind2sub(size(vi),find(vi==num));

asta = zeros(T+1,M);
acon = zeros(T+1,N,M);
c = zeros(T+1,N,M);
v = zeros(T+1,N,M);            % set v0 as a n*m matrix

s=opts;
i=optks;
ks=kmat(optks,1);
theta =(1/(T-s))*1.05;

t=T+1;
  while t>0
        for j=1:M
          if t>T                                                        % find optimal VT
            acon(t,i,j)= amat(j,1);                                     % try grid AT
            c(t,i,j)=((1+r)*acon(t,i,j)-theta*(1-lambda)*ks+y)/p;       % calculate optimal CT
            if c(t,i,j)<0
            v(t,i,j)=-999;
            else
            v(t,i,j)=(1/rho)*(c(t,i,j)^alpha+beta*(ki+ks)^alpha)^(rho/alpha);       % calculate the optimal 
            end
           
           elseif t>s+1                                 % find optimal As+2 to AT
             asta(t,j)= amat(j,1)                       % try grid A
             acon(t,i,j)=optsearch(1,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta);                        % find optimal At+1
             c(t,i,j)=((1+r)*asta(t,j)-acon(t,i,j)-theta*(1-lambda)*ks+y)/p;        % calculate optimal (Cs+1-CT-1)
             v(t,i,j)=utility1(alpha,beta,lambda,rho,p,y,r,ki,amat,T,acon(t,i,j),asta,ks,i,j,t,v,c,theta);            % calculate optimal aggregate utility from t    
 
           elseif t>s                                   % find optimal As+1
              asta(t,j)= amat(j,1)                      % try grid As 
              acon(t,i,j)=optsearch(2,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta);                              % find optimal As+1
              c(t,i,j)=((1+r)*asta(t,j)-acon(t,i,j)-lambda*ks+y)/p                  % calculate optimal Cs
              v(t,i,j)=utility2(alpha,beta,lambda,rho,p,y,r,ki,amat,acon(t,i,j),asta,ks,i,j,t,v);  
  
              elseif t>1                                % find optimal A2 to As; V1
                 asta(t,j)= amat(j,1)                                               % try grid At
                 acon(t,i,j)= optsearch(3,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta)                          % find optimal At+1
                 c(t,i,j)= ((1+r)*asta(t,j)-acon(t,i,j)+y)/p;                       % calculate optimal Ct(C1-Cs-1)
                 v(t,i,j)= utility3(alpha,beta,rho,p,y,r,ki,amat,acon(t,i,j),asta,i,j,t,v);                                  % calculate optimal aggregate utility from t
 
              else
                 acon0= optsearch(4,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta);        % find optimal A1; V0(vi)
                 c0=((1+r)*ai-aconi(s,i)+y)/p;                                                                % find optimal C0   
                 v0= utility4(alpha,beta,rho,p,y,r,ki,ai,amat,aconi(s,i),i,t,v);                  % the optimal value of aggregate utility form t=0
         end
       end
    t=t-1;
 end
   
apath=zeros(T,1);
cpath=zeros(T,1);

apath(1,1)=acon0;
cpath(1,1)=c0;

for t=2:T
alo=max(sum(apath(t-1,1)>amat),1);
ahi=alo+1;
apath(t,1)=acon(t,optks,alo)+(apath(t-1,1)-amat(alo))*(acon(t,optks,ahi)-acon(t,optks,alo))/(amat(ahi)-amat(alo));
cpath(t,1)=c(t,optks,alo)+(apath(t-1,1)-amat(alo))*(c(t,optks,ahi)-c(t,optks,alo))/(amat(ahi)-amat(alo));
end

elapsedTime=tic