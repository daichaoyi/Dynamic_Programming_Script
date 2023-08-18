function tempa=optsearch(stage,alpha,beta,lambda,rho,p,y,r,ki,ai,amat,T,asta,ks,i,j,t,v,c,theta)
 
if stage==1

    grid=0
    tempa=0;
    tempresult=utility1(alpha,beta,lambda,rho,p,y,r,ki,amat,T,0,asta,ks,i,j,t,v,c,theta);

    while grid<14.95
       grid=grid+0.1;
       result=utility1(alpha,beta,lambda,rho,p,y,r,ki,amat,T,grid,asta,ks,i,j,t,v,c,theta);
       if  result>tempresult;
           tempresult=result;
           tempa=grid;
       end
      end


elseif stage==2
    
       grid=0
       tempa=0;
       tempresult=utility2(alpha,beta,lambda,rho,p,y,r,ki,amat,0,asta,ks,i,j,t,v);

       while grid<14.95
        grid=grid+0.1;
         result=utility2(alpha,beta,lambda,rho,p,y,r,ki,amat,grid,asta,ks,i,j,t,v);
         if  result>tempresult;
             tempresult=result;
             tempa=grid;
         end
       end
       
elseif stage==3
    
       grid=0
       tempa=0;
       tempresult=utility3(alpha,beta,rho,p,y,r,ki,amat,0,asta,i,j,t,v);

       while grid<14.95
        grid=grid+0.1;
         result=utility3(alpha,beta,rho,p,y,r,ki,amat,grid,asta,i,j,t,v);
         if  result>tempresult;
             tempresult=result;
             tempa=grid;
         end
       end
       
elseif stage==4
    
       grid=0
       tempa=0;
       tempresult=utility4(alpha,beta,rho,p,y,r,ki,ai,amat,0,i,t,v);

       while grid<14.95
        grid=grid+0.1;
         result=utility4(alpha,beta,rho,p,y,r,ki,ai,amat,grid,i,t,v);
         if  result>tempresult;
             tempresult=result;
             tempa=grid;
         end
       end

 end


        

        

