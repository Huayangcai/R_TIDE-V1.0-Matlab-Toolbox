function yout= Rtide_predict(z,q,t,fu,cof,b,v,f,Qc,fband)
%predict water level by the cofficent b of sin and cos,the optimized cof and
%the driving Q  
%f(:)=1;v(:)=0;
M=length(z);
k=1;
      pl=cof(k,1);TauQ=fix(cof(k,2));
        iq=1:M-TauQ;
        iz=iq+TauQ;
        z1=z(iz,k);m1=find(~isnan(z1));iz1=m1+TauQ;
        time=t(iz,1);m=find(isnan(z1));z1(m)=[];time(m)=[];
        q1=q(iq);q1(m)=[];

     for i=1:length(fband)
       m=find(fu>fband(i,1) & fu<fband(i,2));
       if isempty(m)
           continue
       end
         pr(m)=cof(k,i+2);
       x1(:,m)=[f(m,ones(1,length(time)))]'.*cos(time*2*pi*[fu(m)]'*24+[v(m,ones(1,length(time)))]');
       x2(:,m)=[f(m,ones(1,length(time)))]'.*sin(time*2*pi*[fu(m)]'*24+[v(m,ones(1,length(time)))]');
       Q1=(q1.^cof(k,i+2));
       x3(:,m)=Q1(:,ones(1,length(m))).*x1(:,m);
       x4(:,m)=Q1(:,ones(1,length(m))).*x2(:,m);
     end
     m=find(q1>Qc);x1(m,:)=0;x2(m,:)=0;x3(m,:)=0;x4(m,:)=0;
       xx=[ones(length(x1),1) q1.^pl x1 x2 x3 x4];
       y=z1;%regression model to determine coefficients d and a
       y1=xx*b; %reconstruct water level y1
       si(k)=std(y-real(y1)) %RMSE
     percent(k)=100-sum((y-real(y1)).^2)/sum(y.^2)*100 %R2:Correlation coefficient
     yout(iz1,k)=y1; %predicted water level
%      jj=1:1000;
%      plot(time(jj),z1(jj),time(jj),y1(jj),'r')
end

