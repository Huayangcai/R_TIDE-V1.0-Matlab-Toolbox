function stdmin=Qztidefun(q,z,t,fu,fband,cof,f,v,Qc)
     k=1;
     M=length(q);
     n=length(fu);% number of major tidal constituents 
   %  f(1:n,1:M)=1;v(1:n,1:M)=0;Qc=66600;
    pl=cof(k,1);TauQ=fix(cof(k,2));
        iq=1:M-TauQ;
        iz=iq+TauQ;q1=q(iq);
        z1=z(iz,k);
        st(iz,k)=nan;ft(iz,k)=nan;y2(iz,k)=nan;Eta(iz,1:n,k)=nan;Phi(iz,1:n,k)=nan;
        time=t(iz,1);m=find(isnan(z1));z1(m)=[];time(m)=[];
        q1(m)=[];
        m=find(q1>Qc(1,1));ff(1:length(q1),1)=1;ff(m)=Qc(1,1)./q1(m);%+1.7e-6*(q1(m)-Qc);

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
     m=find(q1>Qc(1,2));x1(m,:)=0;x2(m,:)=0;x3(m,:)=0;x4(m,:)=0;
       xx=[ones(length(x1),1) q1.^pl  x1 x2 x3 x4];
       y=z1;b1 = regress(y,xx);%regression model to determine coefficients d and a
       y1=xx*b1; %reconstruct water level y1
     stdmin=std(y1-z1)
end
    
