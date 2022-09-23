function [y1,y2,yok,yms]= Rtide_discrete(z,q,t,sname)
%discrete diual and semi-diual signal
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
        m=find(q1>Qc);ff(1:length(q1),1)=1;ff(m)=Qc./q1(m)+1.7e-6*(q1(m)-Qc);
        
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
       m=find(fu>fband(2,1) & fu<fband(2,2));
       xx1=[x1(:,m) x2(:,m) x3(:,m) x4(:,m)];
       bb1=[b(m+2);b(n+m+2);b(2*n+m+2);b(3*n+m+2)];
       y1=xx1*bb1;
       m=find(fu>fband(3,1) & fu<fband(3,2));
       xx2=[x1(:,m) x2(:,m) x3(:,m) x4(:,m)];
       bb2=[b(m+2);b(n+m+2);b(2*n+m+2);b(3*n+m+2)];
       y2=xx2*bb2;
       %o1\k1(in this case,sixth and eigth constituent)
       m=[6 8];
       xok=[x1(:,m) x2(:,m) x3(:,m) x4(:,m)];
       bb3=[b(m+2);b(n+m+2);b(2*n+m+2);b(3*n+m+2)];
       yok=xok*bb3;
       %m2\s2(in this case,15/17 constituent)
       m=[15 17];
       xms=[x1(:,m) x2(:,m) x3(:,m) x4(:,m)];
       bb3=[b(m+2);b(n+m+2);b(2*n+m+2);b(3*n+m+2)];
       yms=xms*bb3;
ii=200:500;
plot(t(ii),yms(ii),t(ii),yok(ii))
save con35.mat y1 y2 t yok yms

end