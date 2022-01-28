% R_DEMO - demonstration of capabilities used in Yangtze River Estuary,

% The parameters input:
% xin1 —— water level data used for harmonic analysis
% Q1 —— discharge data used for harmonic analysis
% T1 —— the corresponding time series of the input data in term of 'datenum'
% lat1 —— the latitude of the selected station
% ray —— Rayleigh criteria, the default value is 1, which indicates that the Rayleigh criteria is used to select tidal constituents. Otherwise, Rayleigh criteria is not used.
% synth —— signal noise ratio, the default value is 10. You can set it depending on your own time series
% Qc1 —— the critical discharge beyond which the tide is vanishing, the default value is the maximum of the Q variable.  Generally, it can be set to be the value corresponding with a negligible tidal range last for more than 2 days
% twin —— window spectrum, the default value is 366.
% sname —— the name of the selected station
% ipso —— the method you used for harmonic analysis. If ipso=1, it will invoke PSO to optimize and save the optimized results. If ipso=2, it will directly invoke the optimized file derived last time. So if it is the first time to use it, this parameter should be set as 1.
% ipre —— the method you used for prediction. If ipre=1, it will save the coefficients for prediction and then invoke the program, namely, Rtide_predict.m;. If ipre=2, it will directly predict. Divide the data series into 2 parts, one of which is used to derive coefficient for prediction, another one is used to predict. Users should set the length of data series in the file, namely, Rtide_pre.m. So if it is the first time to use it, this parameter should be set as 1.

% The parameters output:
% nameu —— the name of the selected tidal constituents
% fu —— the frequency of the selected tidal constituents (/h)
% yout —— the reconstructed water level derived from R_TIDE, consisting of st and ft
% st —— the reconstructed water level derived by the residual water level model
% ft —— the reconstructed water level derived by the tidal-fluvial model
% Eta —— the time-dependent amplitude of each tidal constituent derived from R_TIDE
% Phi —— the time-dependent phase of each tidal constituent derived from R_TIDE
% precent —— the correlation coefficient between reconstructed and observed water levels
% si —— the Root Mean Square Error (RMSE) between reconstructed and observed water levels
% cof —— the regression coefficients adopted for each tidal constituent bands derived from R_TIDE

% Reference:
% Contract: Huayang Cai (E-mail: caihy7@mail.sysu.edu.cn); Suying Ou (E-mail: ousuying@mail.sysu.edu.cn)
%%
clc,clear
close all

load Data_Yangtze_river.mat   % Load the example.
mm=find(ZQ(:,1)==datenum(2009,1,1,0,0,0));
ii1=1:mm-1; % select the time series for calibration （2002-2008）
T1=ZQ(ii1,1); % time series from 2002 to 2008
Q1=ZQ(ii1,end); % river discharge from 2002 to 2008
Z1=ZQ(ii1,2:end-1); % water levels from 2002 to 2008
Qc1=max(Q1); % the default value of critical river discharge
ray=1; % Rayleigh criteria
synth=10; % signal-to-noise ratio
twin=366; % days for spectrum window
ipso=2; 
% if ipso=1, using PSO to optimize and save the optimized results.
% if ipso=2, directly using the optimized file derived from PSO last time.
% if ipso=3, using Fmincon function to optimize and save the optimized results.
% if ipso=4, directly using the optimized file derived from Fmincon function last time.
ipre=1;
% If ipre=1, save the coefficients for prediction, invoke Rtide_predict.m;
% If ipre=2, directly predict. Divide the time series into 2 parts, one of which is used to derive coefficient, the one is used for prediction.
% Users should set the length of time series in Rtide_pre.m.

lat=[32.04044;31.9260;32.2293;32.0518;31.7026;31.3291]; % latitude of the tidal gauging stations
% note: different from T_tide and NS_tide, the variables are not conveyed via invoking function, but assigned values directly in the main program
% This program does not provide more options such as the solutions to system of linear algebraic equations, and does not add analysis of inference tide

% end of the input variables

%% calibration
% harmonic analysis
for i=1:size(Z1,2)
    xin1=Z1(:,i); % read the water level of selected time series in each tidal gauging station
    sname=char(stname(i));% read the name of tidal gauging station
    lat1=lat(i); %  read the latitude of tidal gauging station
    dH=0.001;% The critical tidal range, if dH=0.001, the corresponding river discharge is the critical discharge Qc. If the river discharge exceeds this value, there is no tidal signal in this station.    
    Qc1=R_Qtidal(xin1,Q1,T1,dH,Qc1); % calculate the critical river discharge Qc
    [nameu,fu,yout,st,ft,Eta,Phi,percent,si,cof]=R_tide(xin1,Q1,T1,lat1,ray,synth,Qc1,twin,sname,ipso,ipre);% run R_tide program
    if (ipso==1 || ipso==2)
        fname1=['R_TIDE_calibration_for_PSO_' sname '.mat'];
    end
    if (ipso==3 || ipso==4)
        fname1=['R_TIDE_calibration_for_fmin_' sname '.mat'];
    end
    save(fname1,'T1','xin1','Q1','st','ft','percent','Eta','Phi','cof','si','yout','Qc1','nameu'); % save the results   
    
    N=1;
    n=length(fu);% major tidal constituents 
    M1=length(Q1);
    pl=cof(1,1);TauQ=fix(cof(1,2)); 
    iq1=1:M1-TauQ;
    iz1=iq1+TauQ;
    z1=xin1(iz1,1);
    time1=T1(iz1,1);
    QQ1=Q1(iq1);
    for j=1:n    % coefficient matrix
        x11(:,j)=cos(time1*2*pi*fu(j)*24);
        x21(:,j)=sin(time1*2*pi*fu(j)*24);
        x31(:,j)=QQ1.^pl.*x11(:,j);
        x41(:,j)=QQ1.^pl.*x21(:,j);
    end
    xx=[ones(length(x11),1) QQ1.^pl  x11 x21 x31 x41];
    y=z1;
    b1 = regress(y,xx);
    y1=xx*b1; % regression model to determine coefficients d and a
    si=nanstd(y-real(y1)); % RMSE
    percent=nanvar(real(y1))/nanvar(y)*100; % R^2 denotes the Correlation coefficient 
    st(iz1,1)=[ones(length(x11),1) QQ1.^pl ]*b1(1:2); % reconstructed s(t)
    ft(iz1,1)=[x11 x21 x31 x41]*b1(3:end); % reconstructed f(t)
    y2(iz1,1)=y1; % regression to obtain water level

%%
      figure1=figure;
      plot(time1,y,time1,y1,'r',time1,st(iz1,1),'g')
      text(time1(2),max(y1)*0.9,num2str(percent,3))
      dateFormat =2;
      %title(str(k))
      datetick('x',dateFormat)
      xlim([min(time1) max(time1)])
%%
% compute the time-dependent amplitudes and phases of n tidal constituents

  for kk=1:n  %the number of tidal constituents
      Ak=sqrt(b1(kk+2).^2+b1(n+kk+2).^2); Bk1=sqrt(b1(2*n+kk+2).^2+b1(3*n+kk+2).^2);
      Bk=QQ1.^pl.*sqrt(b1(2*n+kk+2).^2+b1(3*n+kk+2).^2);
      Alphak=atan2(b1(n+kk+2),b1(kk+2));
      Betak=atan2(b1(3*n+kk+2),b1(2*n+kk+2));
      Zk1=0.5*(Ak*cos(Alphak)+Bk.*cos(Betak));
      Zk2=0.5*(Ak*sin(Alphak)+Bk.*sin(Betak));
      Zfk=Zk1+Zk2*j;
      Zk=Zk1-j*Zk2;
      Eta2(iz1,kk)=abs(Zk)+abs(Zfk);
      Phi2(iz1,kk)=90-atan2d(Zk2,Zk1);
      m=find(Phi2(iz1,kk)<0);Phi2(m,kk)=Phi2(m,kk)+360;
      ab(1:7,kk)=[Ak Bk1 Alphak Betak b1(1) b1(2) b1(3)];

  end
jj=0;
   for k=1:24:size(xin1,1)-24-floor(cof(2))
    k1=k; jj=jj+1;   k2=k1+23; kk1=k1:k2;
    zm(jj,1)=mean(xin1(kk1,1));
    ym(jj,1)=mean(yout(kk1,1));
    stm(jj,1)=mean(st(kk1,1));
    ftm(jj,1)=mean(ft(kk1,1));
    qu(jj,1)=mean(QQ1(kk1));
    for kk=1:n
    Eta(jj,kk)=mean(Eta2(kk1,kk));
    Phi(jj,kk)=mean(Phi2(kk1,kk));
    end
    yy1(iz1,1)=y1;
   end
 % compute variances for s(t)_zr and f(t)_ztr through hourly data and daily averaged data
  varhour(k,1)=var(y); 
  varhour(k,2)=var(st(:,1)); 
  varhour(k,3)=var(ft(:,1)); 
  varday(k,1)=var(zm(:,1)); 
  varday(k,2)=var(stm(:,1)); 
  varday(k,3)=var(ftm(:,1)); 
  clear z1 u1 timet1 x11 x21 x31 x41 x51 x61
    if (ipso==1 || ipso==2)
        fname2=['R_TIDE_outputforper_for_PSO_' sname '.mat'];
    end
    if (ipso==3 || ipso==4)
        fname2=['R_TIDE_outputforper_for_fmin_' sname '.mat'];
    end
    save (fname2,'xin1','y2','Eta','Phi','ab','cof','b1','percent','st','ft','QQ1')
    %% validation
    % In this part, we use the coefficients derived from calibration period
    % as well as the river discharge records in the rest time series to
    % reconstruct the water level series, which reveal the influence of the
    % alteration in river discharge, ignoring the geometric change
 % ***************************************************************************************** 
    
 ii2=mm:size(ZQ,1); % time series from 2009 to 2012
 z2=ZQ(ii2,2:end-1);
    xin2=z2(:,i);% water level from 2009 to 2012
    T2=ZQ(ii2,1);% time series from 2009 to 2012
    Q2=ZQ(ii2,end); % river discharge record from 2009 to 2012
    M2=length(Q2);
    
    pl=cof(1,1);TauQ=fix(cof(1,2)); 
    iq2=1:M2-TauQ;
    iz2=iq2+TauQ;
    z2=xin2(iz2,1);
    time2=T2(iz2,1);
    QQ2=Q2(iq2);
    for j=1:n
         x12(:,j)=cos(time2*2*pi*fu(j)*24);
         x22(:,j)=sin(time2*2*pi*fu(j)*24);
         x32(:,j)=QQ2.^pl.*x12(:,j);
         x42(:,j)=QQ2.^pl.*x22(:,j);
    end
     clear j
       xx=[ones(length(x12),1) QQ2.^pl  x12 x22 x32 x42];
       yy2=xx*b1; %regression model to determine coefficients d and a
       Bias=nanmean(z2-yy2);
       si=nanstd(z2-real(yy2)); %RMSE
       %st(iz,k)=[ones(length(x1),1) q1.^pl ]*b1(1:2); %reconstructed s(t)
       %ft(iz,k)=[x1 x2 x3 x4]*b1(3:end); %reconstructed f(t)
       YY2(iz2,1)=yy2; %regression to obtain water level
      percent=nanvar(real(yy2))/nanvar(z2)*100; %R2:Correlation coefficient 
     for kk=1:n  
      Ak=ab(1,kk);
      Bk=QQ2.^pl*ab(2,kk);
      Alphak=ab(3,kk);
      Betak=ab(4,kk);
      Zk1=0.5*(Ak*cos(Alphak)+Bk.*cos(Betak));
      Zk2=0.5*(Ak*sin(Alphak)+Bk.*sin(Betak));
      Zfk=Zk1+Zk2*j;
      Zk=Zk1-j*Zk2;
      Etay(iz2,kk)=abs(Zk)+abs(Zfk);
      Phiy(iz2,kk)=90-atan2d(Zk2,Zk1);
      %m=find(Phi2(iz,kk)<0);Phi2(m,kk)=Phi2(m,kk)+360;
  end
      bb1(1:3)=ab(5:7,1,1);
      sty(iz2,1)=[ones(length(z2),1) QQ2.^pl ]*bb1(1:2)';
      jj=0;
      % calculate the dailymean value
      m=find(Etay==0);Etay(m)=nan;Etam=reshape(nanmean(reshape(Etay(1:end-20,:),24,[])),[],length(fu));
      m=find(Phiy==0);Phiy(m)=nan;Phim=reshape(nanmean(reshape(Phiy(1:end-20,:),24,[])),[],length(fu));
      m=find(xin2==0);xin2(m)=nan;ym=reshape(nanmean(reshape(xin2(1:end-20,:),24,[])),[],1);
      m=find(YY2==0);YY2(m)=nan;yym=reshape(nanmean(reshape(YY2(1:end-20,:),24,[])),[],1);
      if (ipso==1 || ipso==2)
        fname3=['R_TIDE_validation_for_PSO_' sname '.mat'];
    end
    if (ipso==3 || ipso==4)
        fname3=['R_TIDE_validation_for_fmin_' sname '.mat'];
    end
    save(fname3,'Etay','Etam','Phiy','Phim','xin2','ym','yym','nameu')
    clear z2 u2 timet1 x12 x22 x32 x42 x52 x62
end
