function [percent,si,percent1,si1]= Rtide_pre(z,q,t,fu,cof,v,f,Qc,fband)
%non-stationary harmonic analysis by the optimized cof and the driving Q
%f(:)=1;v(:)=0;

n=length(fu);% number of major tidal constituents 
[M,N]=size(z); %N,number of stations;M,length of z and q data;
for k=1:N 
      pl=cof(k,1);TauQ=fix(cof(k,2));
        iq=1:M-TauQ;
        iz=iq+TauQ;
        z1=z(iz,k);m1=find(~isnan(z1));iz1=m1+TauQ;
        st(iz,k)=nan;ft(iz,k)=nan;y2(iz,k)=nan;Eta(iz,1:n,k)=nan;Phi(iz,1:n,k)=nan;
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
      ii=1:61320;
       xx=[ones(length(x1(ii,:)),1) q1(ii,:).^pl  x1(ii,:) x2(ii,:) x3(ii,:) x4(ii,:)];
       y=z1(ii);b1 = regress(y,xx);y1=xx*b1; %regression model to determine coefficients d and a
       si(k)=std(y-real(y1)); %RMSE
     percent(k)=(1-var(y1-z1(ii))/var(z1(ii)))*100; %R2:Correlation coefficient 
%predict
      ii1=61321:length(q1);
       xx1=[ones(length(x1(ii1,:)),1) q1(ii1,:).^pl  x1(ii1,:) x2(ii1,:) x3(ii1,:) x4(ii1,:)];
       y2=xx1*b1;
       si1(k)=std(z1(ii1)-y2) %RMSE
     percent1(k)=(1-var(y2-z1(ii1))/var(z1(ii1)))*100 %R2:Correlation coefficient 


    clear z1 u1 timet1 x1 x2 x3 x4 x5 x6 m

end