function u1=utility1(alpha,beta,lambda,rho,p,y,r,ki,amat,T,a,asta,ks,i,j,t,v,c,theta)

c1(t,i,j)=((1+r)*asta(t,j)-a-theta*(1-lambda)*ks+y)/p;
alo=max(sum(a>amat),1);
ahi=alo+1;

if  c1(t,i,j)<0 || c(T+1,i,alo)<0
    u1=-999-9*abs(c1(t,i,j));
else 
    gg= v(t+1,i,alo)+(a-amat(alo))*(v(t+1,i,ahi)-v(t+1,i,alo))/(amat(ahi)-amat(alo));
    u1=(1/rho)*(c1(t,i,j)^alpha+beta*(ki+ks)^alpha)^(rho/alpha)+beta*gg;
end    

