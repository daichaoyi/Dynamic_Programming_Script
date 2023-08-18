function u2=utility2(alpha,beta,lambda,rho,p,y,r,ki,amat,a,asta,ks,i,j,t,v)     

c2(t,i,j)=((1+r)*asta(t,j)-a-lambda*ks+y)/p;   
alo=max(sum(a>amat),1);
ahi=alo+1;

if  c2(t,i,j)<0
    u2=-999-9*abs(c2(t,i,j));
else
    gg=v(t+1,i,alo)+(a-amat(alo))*(v(t+1,i,ahi)-v(t+1,i,alo))/(amat(ahi)-amat(alo));
    u2=(1/rho)*(c2(t,i,j)^alpha+beta*(ki+ks)^alpha)^(rho/alpha)+beta*gg;
end
 