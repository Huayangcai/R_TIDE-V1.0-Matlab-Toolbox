function [z1,y2,percent1,si1,st,ft]= Rtide_pre(z,q,t,sname)
%non-stationary harmonic analysis by the optimized cof and the driving Q
%f(:)=1;v(:)=0;
  fname1=['harmcofforpre_' sname '_2002_2008' '.mat'];
  load(fname1)
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
        m=find(q1>Qc(k,1));
        ff(1:length(q1),1)=1;
        ff(m,1)=Qc(k,1)./q1(m,1);%+1.7e-6*(q1(m)-Qc);
        
     for i=1:length(fband)
       m=find(fu>fband(i,1) & fu<fband(i,2));
       if isempty(m)
           continue
       end
       x1(:,m)=[f(m,ones(1,length(time)))]'.*cos(time*2*pi*[fu(m)]'*24+[vu(m,ones(1,length(time)))]');
       x2(:,m)=[f(m,ones(1,length(time)))]'.*sin(time*2*pi*[fu(m)]'*24+[vu(m,ones(1,length(time)))]');
       Q1=((q1.*ff).^cof(k,i+2));
       x3(:,m)=Q1(:,ones(1,length(m))).*x1(:,m);
       x4(:,m)=Q1(:,ones(1,length(m))).*x2(:,m);
     end
       xx=[ones(length(x1),1) q1.^pl  x1 x2 x3 x4];
       y2=xx*b;
      si1=nanstd(z1-real(y2)); % RMSE
      percent1=100-nanvar(z1-real(y2))/nanvar(z1)*100; % R^2 denotes the Correlation coefficient 
    st(iz1,1)=[ones(length(x1),1) q1.^pl ]*b(1:2); % reconstructed s(t)
    ft(iz1,1)=[x1 x2 x3 x4]*b(3:end); % reconstructed f(t)
%%
      figure1=figure;
      y=z1;
      plot(time,y,time,y2,'r',time,st(iz1,1),'g')
      title(['R^2= ' num2str(percent1/100,3)],'fontname','Times New Roman')
      xlabel('Date')
      ylabel('\itZ \rm(m)')
      dateFormat =2;
      %title(str(k))
      datetick('x',dateFormat)
      xlim([min(time) max(time)])
      set(gca,'fontname','Times New Roman')
      legend('\itZ_{\rmObs.}','\itZ_{\rmSim.}','s(t)')

 %   clear z1 u1 timet1 x1 x2 x3 x4 x5 x6 m

end
