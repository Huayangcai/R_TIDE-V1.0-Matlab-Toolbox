      function fx=Rtide_sixcof(X,q,z,t,fu,fband,id,v,f,Qc)
   %least square method for R2 or RMSE
     [N,n]=size(X);M=length(q);
     for k=1:N
      pl=X(k,1);TauQ=round(X(k,2));  pr(1:n-2)=X(k,3:end);
          iq=1:M-TauQ;
          iz=iq+TauQ;
      z1=z(iz,1);
      time=t(iz,1);m=find(isnan(z1));z1(m)=[];time(m)=[];
     q1=q(iq);q1(m)=[];
       ff(1:length(q1),1)=1; m=find(q1>Qc(1,1));
        ff(m,1)=Qc(1,1)./q1(m);
     %find low freq constituent
     for i=1:length(fband)
       m=find(fu>fband(i,1) & fu<fband(i,2));
       x1(:,m)=[f(m,ones(1,length(time)))]'.*cos(time*2*pi*[fu(m)]'*24+[v(m,ones(1,length(time)))]');
       x2(:,m)=[f(m,ones(1,length(time)))]'.*sin(time*2*pi*[fu(m)]'*24+[v(m,ones(1,length(time)))]');
       Q1=((q1.*ff).^pr(i));
       x3(:,m)=Q1(:,ones(1,length(m))).*x1(:,m);
       x4(:,m)=Q1(:,ones(1,length(m))).*x2(:,m);
       pr1(m)=pr(i);
     end
     m=find(q1>Qc(1,2));x1(m,:)=0;x2(m,:)=0;x3(m,:)=0;x4(m,:)=0; %%(杩欓噷鎴戞敼浜嗭紙1锛?2锛夛紝---leebok)
       xx=[ones(length(x1),1) q1.^pl  x1 x2 x3 x4];
       y=z1;b = regress(y,xx);y1=xx*b; %regression model to determine coefficients d and a

           si2(k)=std(y-real(y1)); %RMSE
           r2(k)=var(real(y1)-mean(y))/var(y)*100; %R2:Correlation coefficient 

         clear z1 time q1 t1 x1 x2 x3 x4 x5 x6
     end
     
              
       if id==1
           fx=r2;
       else
           fx=si2;
       end
