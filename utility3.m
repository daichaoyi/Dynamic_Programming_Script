function  u3=utility3(alpha,beta,rho,p,y,r,ki,amat,a,asta,i,j,t,v)

c3(t,i,j)=((1+r)*asta(t,j)-a+y)/p;  
    alo=max(sum(a>amat),1);
    ahi=alo+1;
if c3(t,i,j)<0
   u3=-999999-9*abs(c3(t,i,j));
else
    gg=v(t+1,i,alo)+(a-amat(alo))*(v(t+1,i,ahi)-v(t+1,i,alo))/(amat(ahi)-amat(alo));
    u3=(1/rho)*(c3(t,i,j)^alpha+beta*(ki)^alpha)^(rho/alpha)+beta*gg;
end
