function [st,ft,yout,Eta,Phi,percent,si,b]= Rtide_harmonic(z,q,t,fu,cof,v,f,Qc,fband)
%non-stationary harmonic analysis by the optimized cof and the driving Q

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
        m=find(q1>Qc(1,1));ff(1:length(q1),1)=1;ff(m)=Qc(1,1)./q1(m);

     for i=1:length(fband)
       m=find(fu>fband(i,1) & fu<fband(i,2));
       if isempty(m)
           continue
       end
         pr(m)=cof(k,i+2);
       x1(:,m)=[f(m,ones(1,length(time)))]'.*cos(time*2*pi*[fu(m)]'*24+[v(m,ones(1,length(time)))]');
       x2(:,m)=[f(m,ones(1,length(time)))]'.*sin(time*2*pi*[fu(m)]'*24+[v(m,ones(1,length(time)))]');
       Q1=((q1.*ff).^cof(k,i+2));
       x3(:,m)=Q1(:,ones(1,length(m))).*x1(:,m);
       x4(:,m)=Q1(:,ones(1,length(m))).*x2(:,m);
     end
     m=find(q1>Qc(1,2));x1(m,:)=0;x2(m,:)=0;x3(m,:)=0;x4(m,:)=0;  %%(è¿™é‡Œæ”¹äº†ï¼?1ï¼?2ï¼‰â?”â?”â?”â?”â?”â?”leebok)
       xx=[ones(length(x1),1) q1.^pl  x1 x2 x3 x4];
       y=z1;b1 = regress(y,xx);%regression model to determine coefficients d and a
       y1=xx*b1; %reconstruct water level y1
       si(k)=std(y-real(y1)); %RMSE
%     percent(k)=100-sum((y-real(y1)).^2)/sum(y.^2)*100; %R2:Correlation coefficient 
     percent(k)=100-var((y-real(y1)))/var(y)*100; %R2:Correlation coefficient 
     st(iz1,k)=[ones(length(x1),1) q1.^pl ]*b1(1:2); %reconstructed s(t)
     ft(iz1,k)=[x1 x2 x3 x4]*b1(3:end); %reconstructed f(t)
     yout(iz1,k)=y1; %regression water level
     b(:,k)=b1;
%%
     % plot test figure;
%       plot(time,y,time,y1,'r',time,st(iz1,k),'g')
%       text(time(2),max(y1)*0.9,num2str(percent(k),3))
%       dateFormat =2;
%       %title(str(k))
%       datetick('x',dateFormat)
%       xlim([min(time) max(time)])
%%
% compute the time-dependent amplitudes and phases of n tidal constituents
  for kk=1:n  % n constituent
      Ak=sqrt(b1(kk+2).^2+b1(n+kk+2).^2);
      Bk=(q1.*ff).^pr(kk).*sqrt(b1(2*n+kk+2).^2+b1(3*n+kk+2).^2);
      Alphak=atan2(b1(n+kk+2),b1(kk+2));
      Betak=atan2(b1(3*n+kk+2),b1(2*n+kk+2));
      Zk1=0.5*(Ak*cos(Alphak)+Bk.*cos(Betak));
      Zk2=0.5*(Ak*sin(Alphak)+Bk.*sin(Betak));
      Zfk=Zk1+Zk2*j;
      Zk=Zk1-j*Zk2;
      Eta(iz1,kk,k)=abs(Zk)+abs(Zfk);
      Phi(iz1,kk,k)=atan2d(Zk2,Zk1);
%       [minx,ix]=min(Eta(iz1,kk,k));
 %      if ix>iq/3
%            m=find(q1>Qc);
%         Eta(iz1(m),kk,k)=NaN;
%         Phi(iz1(m),kk,k)=NaN;
 %      end
%     Eta(m,kk,k)=sqrt((b1(kk+2)+q1(m).^pr(kk).*b1(2*n+kk+2)).^2+(b1(n+kk+2)+q1(m).^pr(kk).*b1(3*n+kk+2)).^2);
%     Phi(m,kk,k)=atan2d((b1(n+kk+2)+q1(m).^pr(kk).*b1(3*n+kk+2)),(b1(kk+2)+q1(m).^pr(kk).*b1(2*n+kk+2)));
  end

    clear z1 u1 timet1 x1 x2 x3 x4 x5 x6 m

end