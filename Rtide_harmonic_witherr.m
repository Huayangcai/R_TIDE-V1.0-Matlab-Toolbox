function [st,ft,yout,yout_snr,percent,si,b,Eta,Phi,tidecon]= Rtide_harmonic_witherr(z,q,t,fu,cof,v,f,Qc,fband,synth)
%non-stationary harmonic analysis by the optimized cof and the driving Q

n=length(fu);% number of major tidal constituents 
[M,N]=size(z(:,1)); %N,number of stations;M,length of z and q data;
f(1:n,1:M)=1;v(1:n,1:M)=0;

for k=1:N 
      pl=cof(k,1);TauQ=fix(cof(k,2));
        iq=1:M-TauQ;
        iz=iq+TauQ;
        z1=z(iz,k);m1=find(~isnan(z1));iz1=m1+TauQ;
        st(iz,k)=nan;ft(iz,k)=nan;y2(iz,k)=nan;Eta(iz,1:n,k)=nan;Phi(iz,1:n,k)=nan;
        time=t(iz,1);m=find(isnan(z1));z1(m)=[];time(m)=[];
        q1=q(iq);q1(m)=[];
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
      xx1=[q1.^pl  x1 x2 x3 x4];
      [b2,stats]=robustfit(xx1,y);
      y1=xx*b1; %reconstruct water level y1
       si(k)=std(y-real(y1)); %RMSE
     percent(k)=100-sum((y-real(y1)).^2)/sum(y.^2)*100; %R2:Correlation coefficient 
     st(iz1,k)=[ones(length(x1),1) q1.^pl ]*b1(1:2); %reconstructed s(t)
     ft(iz1,k)=[x1 x2 x3 x4]*b1(3:end); %reconstructed f(t)
     yout(iz1,k)=y1; %regression water level
     b(:,k)=b1;
end
%%
% error estimated,use correlated noise error estimated from NS_tide
cases=300;
    TC1 =xx;
    W = stats.w; %weigths
    xres=y-y1;
    xres_ci = W(1:length(xres)).*xres;
    A = W(:,ones(1,size(TC1,2))).*TC1;
    dum = sum(A,2)+z1; 
    Agd = ~isnan(dum); %locate the non-NaNs
    Pf = fft(A(Agd,:));
    Perr = fft(xres_ci(Agd));
%     Perr = fft(xres(Agd));
    %Daniell avg the spectrum (boxcar avg)
    B = ones(16,1) ;
    B = B./sum(B) ;
    msk = ones(length(Perr),1) ;
%     Perr = conv(abs(Perr).^2,B,'same')./conv(msk,B,'same') ; %for newer Matlab versions
    Perr = conv2(abs(Perr).^2,B,'same')./conv2(msk,B,'same') ; %for older Matlab versions
    covmat = Pf'*((1./Perr(:,ones(1,size(Pf,2)))).*Pf);
    covmat = (covmat + covmat')/2; %symmetrizing the matrix
%     SIGMA = inv(Pf'*(diag(1./Perr)*Pf));
    SIGMA = real(inv(covmat));
    MU = b1';
    if length(SIGMA)~=length(MU)
        warning('Length of SIGMA is different than length of MU');
        for i=1:length(MU)-length(SIGMA)
            SIGMA(end+1,end+1) = 0;
        end
    end
    coefr = mvnrnd(MU,SIGMA,cases-1);
    coefr = [b1'; coefr]; %first row is the calculated coefficients from robustfit

%Coefficients
c=[];         %stage coefficients
c(1,:) = (coefr(:,1))';
ic = 2;       %iterator on the coefficients
id = 1;       %iterator on the d matrix
ii=1:cases;
c(2:2,:) = (coefr(:,2:2));
for i=1:n
    dc(i,:)=coefr(:,i+2);
    dc(i+n,:)=coefr(:,2*n+i+2); %cos
    ds(i,:)=coefr(:,n+i+2);
    ds(i+n,:)=coefr(:,3*n+i+2); %sin
end
[mc nc] = size(c);
[md nd] = size(dc);
%a+ and a- amplitudes
ap=(dc-1i*ds)/2;
am=(dc+1i*ds)/2;
  epsp=angle(ap)*180/pi;
  epsm=angle(am)*180/pi;
  ap=abs(ap);
  am=abs(am);
aap=ap./[f(:,ones(1,cases)); f(:,ones(1,cases))];	% Apply nodal corrections and
aam=am./[f(:,ones(1,cases)); f(:,ones(1,cases))];	% compute ellipse parameters.

fmaj=aap+aam;                   % major axis
fmin=aap-aam;                   % minor axis

gp=mod( [v(:,ones(1,cases)); v(:,ones(1,cases))]-epsp ,360); % pos. Greenwich phase in deg.
gm=mod( [v(:,ones(1,cases)); v(:,ones(1,cases))]+epsm ,360); % neg. Greenwich phase in deg.

finc= (epsp+epsm)/2;
finc(:,1)=mod( finc(:,1),180 ); % Ellipse inclination in degrees
	
finc=cluster(finc,180); 	% Cluster angles around the 'true' 
                                % angle to avoid 360 degree wraps.

pha=mod( gp+finc ,360); 	% Greenwich phase in degrees.

pha=cluster(pha,360);		% Cluster angles around the 'true' angle
				% to avoid 360 degree wraps.
%Global test
    DpGL=[]; fmajGL=[]; phaGL=[];  %model parameters
    CTS=[]; fmajTS=[]; phaTS=[];   %final time series
    emaj=[]; epha=[]; eCTS=[];     %error bounds
    inn=length(z1);fmajTS(1:n,1:inn)=nan;phaTS(1:n,1:inn)=nan;
 %stage
     CTS = zeros(1,inn);  %sum of all terms, 
     eCTS= zeros(1,inn);  %with corresponding error 'eCTS'
       %1)construct the matrix and generate 95% CI
     CTS = CTS + c(1,1)*ones(1,inn);
     eCTS = eCTS + median(abs(ones(inn,1)*c(1,:)-median(ones(inn,1)*c(1,:),2)*ones(1,cases)),2)'/.6745*1.96;
     CTS = CTS + c(2,1)*(q1.^pl)';
     eCTS = eCTS + median(abs(q1(:,1).^pl*c(2,:)-median(q1(:,1).^pl*c(2,:),2)*ones(1,cases)),2)'/.6745*1.96;
% harmonic
     [md nd] = size(fmaj);
            %1)calculate model parameter and generate 95% CI\
   for i=1:n
       
                DpGL = zeros(inn,cases);
       DpGL = DpGL + ones(inn,1)*((fmaj(i,:)/2.0).*exp(-1i*pha(i,:)*pi/180));
%       DpGL = DpGL + (min(1,(Qc(ij,1)./sum(Q,2)).^Qc(ij,2)))*(fmaj{id,1}(i,:)/2.*exp(-1i*pha{id,1}(i,:)*pi/180));
       DpGL = DpGL + ((q1.*ff).^pr(i)*(fmaj(i+n,:)/2.*exp(-1i*pha(i+n,:)*pi/180)));
       fmajGL = abs(DpGL) + abs(conj(DpGL));
       phaGL  = mod(imag(log(conj(DpGL)))*180/pi, 360);
         fmajTS(i,iz1) = fmajGL(:,1);
          phaTS(i,iz1)  = phaGL(:,1);
%         fmajTS(i,iz1) = mean(fmajGL');
%          phaTS(i,iz1)  = mean(phaGL');
         emaj(i,iz1) = median(abs(fmajGL-median(fmajGL,2)*ones(1,cases)),2)/.6745*1.96;
         epha(i,iz1) = median(abs( phaGL-median( phaGL,2)*ones(1,cases)),2)/.6745*1.96;
         snr(i,iz1) = (fmajTS(i,iz1)./emaj(i,iz1)).^2;
 % snr>synth 
        snrm(i)=mean(snr(i,iz1));
   end
   tidecon.fmaj=fmajTS; tidecon.phi=phaTS; tidecon.emaj=emaj; tidecon.ephi=epha; tidecon.snr=snr;
%%
% % compute the time-dependent amplitudes and phases of n significant tidal constituents
   m=find(snrm>synth)
   xxm=[ones(length(x1),1) q1.^pl  x1(:,m) x2(:,m) x3(:,m) x4(:,m)];
   b3 = regress(z1,xxm);
   yout_s=xxm*b3;
   yout_snr(iz1,1)=yout_s;
   si_snr=std(z1-yout_s)
   n=length(m);
  for kk=1:n  % n constituent
      Ak=sqrt(b3(kk+2).^2+b3(n+kk+2).^2);
      Bk=q1.^pr(m(kk)).*sqrt(b3(2*n+kk+2).^2+b3(3*n+kk+2).^2);
      Alphak=atan2(b3(n+kk+2),b3(kk+2));
      Betak=atan2(b3(3*n+kk+2),b3(2*n+kk+2));
      Zk1=0.5*(Ak*cos(Alphak)+Bk.*cos(Betak));
      Zk2=0.5*(Ak*sin(Alphak)+Bk.*sin(Betak));
      Zfk=Zk1+Zk2*j;
      Zk=Zk1-j*Zk2;
      Eta(iz1,m(kk))=abs(Zk)+abs(Zfk);
      Phi(iz1,m(kk))=atan2d(Zk2,Zk1);
  end


