function  u4=utility4(alpha,beta,rho,p,y,r,ki,ai,amat,a,i,t,v)

c4(t,i)=((1+r)*ai-a+y)/p;
if c4(t,i)<0
   u4=-999999-9*abs(c4(t,i));
else
    alo=max(sum(a>amat),1);
    ahi=alo+1;
    gg=v(t+1,i,alo)+(a-amat(alo))*(v(t+1,i,ahi)-v(t+1,i,alo))/(amat(ahi)-amat(alo));
    u4=(1/rho)*(c4(t,i)^alpha+beta*(ki)^alpha)^(rho/alpha)+beta*gg;
end
