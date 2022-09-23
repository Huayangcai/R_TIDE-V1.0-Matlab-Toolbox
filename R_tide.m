function [nameu,fu,yout,st,ft,Eta,Phi,percent,si,cof,tidecon]=R_tide(xin,Q,t,lat,ray,synth,Qc,twin,sname,IPso,ipre,tsnr,tau);
% R_TIDE Harmonic analysis of a tidal time series modified by non-stationary river discharge Q 
 
dt=1;%interval time
fid=1;
Isnr=0;
stime=t(1); %start time
corr_fs=[0 1e6];
corr_fac=[1  1];
secular='mean';
shallownames=[];
constitnames=[];
centraltime=[];
errcalc='cboot'; %error estimate method
lsq='direct'; %least square method
       % Define inference parameters.
       inf.iname=['P1';'K2'];
       inf.irefname=['K1';'S2'];
       inf.amprat=[.33093;.27215];
       inf.ph=[-7.07;-22.40];
  
%input variable end
 %%
[inn,inm]=size(xin);
if ~(inn==1 | inm==1), error('Input time series is not a vector'); end;%input tidal signal must be column

xin=xin(:); % makes xin a column vector
nobs=length(xin);

if nobs*dt> 18.6*365.25*24,  % Long time series
  longseries=1; ltype='full';
else
  longseries=0; ltype='nodal';
end;
        		
nobsu=nobs-rem(nobs-1,2);% makes series odd to give a center point

%t=stime+dt*[1:nobs]'/24.0;  % Time vector for entire time series,Julian time(datenum time)
%centraltime=stime+dt*[1:nobs]'/24.0/2;
centraltime=t(fix(nobsu/2));
    fband =[0 .03;
          .03 .06;
          .06 .10;
          .10 .14;
          .14 .18;
          .18 .22;
          .22 .26;
          .26 .29;
          .30 .50];
% -------Get the frequencies to use in the harmonic analysis-----------
%combined with Mattes(2013),consider that effect of River discharge on the
%constituents
minres=ray/(dt*nobsu);
% number of tidal constituents without river discharge
[nameu,fu,ju,namei,fi,jinf,jref]=constituents(minres,constitnames,...
                                           shallownames,[],[],centraltime);
%[nameu,fu,ju,namei,fi,jinf,jref]=constituents(minres,constitnames,...
%                                           shallownames,inf.iname,inf.irefname,centraltime);
expcof=1.5;% a default parameters for tidal constituent select
crit=0.1; % eta cofficient,this value can adjust from 0.001~0.5 to define delta sigma,to select constituent mathod refers to Ns_tide
timestep = round((t(end)-t(1))/length(t)*24*60)/60; %in hours
[inn,inm] = size(xin);
if inn<twin*24/timestep   
    win = [];  
else
    win = twin*24/timestep;  
end
nameu1=[]; fu1=[]; ju1=[];

 for jm=1:length(fband) %loop on the frequency intervals (fm>=1 when model is 'fluvial')

        %initial frequencies
   fui = fu(fu>fband(jm,1) & fu<fband(jm,2));
   minres1 = R_rayleigh(t,Q,expcof,fui,win,crit,minres)
  [nameu_ray,fu_ray,ju_ray,namei,fi,jinf,jref]=constituents(minres1,constitnames,...
                                           shallownames,[],[],centraltime);
    i1=fu_ray>fband(jm,1) & fu_ray<fband(jm,2); 
     nameu1 = [nameu1; nameu_ray(i1,:)];
    fu1    = [fu1; fu_ray(i1)];
    ju1    = [ju1; ju_ray(i1)];
 end
  clear nameu fu ju
  nameu=nameu1;fu=fu1;ju=ju1;
mu=length(fu); % # base frequencies
mi=length(fi); % # inferred

% Find the good data points (here I assume that in a complex time 
% series, if u is bad, so is v).

gd=find(isfinite(xin(1:nobsu)));
ngood=length(gd);
fprintf('   Points used: %d of %d\n',ngood,nobs)

%---------------Nodal Corrections-------------------------------------- 						   
% Generate nodal corrections and calculate phase relative to Greenwich. 						   
% Note that this is a slightly weird way to do the nodal corrections,							   
% but is 'traditional'.  The "right" way would be to change the basis							   
% functions used in the least-squares fit above.									   
%note:it's simple,no change
  % Get nodal corrections at midpoint time.	
  [v,u,f]=t_vuf(ltype,(t(end)+t(1))/2,[ju;jinf],lat);									   
  vu=(v+u);%*360; % total phase correction 									   

%%
% vu1=u(:,1:length(ju));f1=f(:,1:length(ju));
%     tc=[ones(length(t),1) f1*[cos((2*pi)*t*fu'*24-vu1)] f1.*sin((2*pi)*t*fu'*24-vu1) ];
     tc=[ones(length(t),1) [cos((2*pi)*t*fu'*24)] sin((2*pi)*t*fu'*24) ];

   %coef=tc(gd,:)\xin(gd);
   coef=regress(xin(gd),tc(gd,:));
% 
  z0=coef(1);
  ap=(coef(2:(1+mu))-i*coef((2+mu):(1+2*mu)))/2;  % a+ amplitudes
  am=(coef(2:(1+mu))+i*coef((2+mu):(1+2*mu)))/2;  % a- amplitudes
% 

% --------------Error Bar Calculations and 95% confidence---------------------------------
  xout=tc*coef;  % This is the time series synthesized from the analysis
  xres=xin-xout; % and the residuals!
%   plot(t,xin,t,xout)
%   sum(xres.^2)/sum((xin-mean(xin)).^2)
%   sum((xout-mean(xin)).^2)/sum((xin-mean(xin)).^2)
   %r2=var(xout-mean(xin))/var(xin) 
%-----select significant tide---------
[snr,tidecon]=error_snr(xres,fu,dt,errcalc,ap,am,f,vu);

%according to SNR  to select significant constituent and decrease number of
%tide constituent
if Isnr==1
      for k=1:length(fband)
          m=find(fu>=fband(k,1)&fu<fband(k,2));
          fum(m)=nan;
          m1=find(snr(m)>synth);%select typical constituents which snr>sy
          if ~isempty(m1)
             xm=max(tidecon(m,1));i2=find(tidecon(m,1)/xm>0.1&snr(m)>synth);
              fum(m(i2))=fu(m(i2));
          else
             [~,m0] =max(tidecon(m,1));
             fum(m(m0))=fu(m(m0));
          end
      end
      m=find(isnan(fum));
fu(m)=[];nameu(m,:)=[];ju(m)=[];
end

%%
%add astro constituent SA
%fu=[fu;0.0001141];nameu=[nameu;'SA  '];f=[f;1];vu=[vu;0];
%PART III£ºcoefficient optimization , using the fmincon in Matlab
    fband =[0 .03;
          .03 .06;
          .06 .10;
          .10 .14;
          .14 .18;
          .18 .22;
          .22 .26;
          .26 .29;
          .30 .50];
      time=datevec(t);
%IPso=2;  %if first run, it must run function to optimize cof;and IPso=2, it have cof,and give harmonic quickly
if IPso==3
  [cof,stmin]=fmin_Q_tide(xin,Q,t,fu,fband,f,vu,Qc,tau); %it would some time to optimize,per is R2
  cof=cof';
  if time(1,1)==time(end,1)
      optfilename=['opti_fmin_' sname '_' num2str(time(1,1)) '.mat'];
  else
      optfilename=['opti_fmin_' sname '_' num2str(time(1,1)) '_' num2str(time(end,1)) '.mat'];
  end
  save(optfilename,'cof','stmin','fband','f')   %so if opti ok,save it,then next test again,use load to instead of run again
end
if IPso==1
    id=2;
  [cof,per,fband]=Rtide_pso(xin,Q,t,fu,id,vu,f,Qc,fband,tau); %it would some time to optimize,per is R2
 if time(1,1)==time(end,1)
      optfilename=['opti_pso_' sname '_' num2str(time(1,1)) '.mat'];
  else
      optfilename=['opti_pso_' sname '_' num2str(time(1,1)) '_' num2str(time(end,1)) '.mat'];
  end
  save(optfilename,'cof','per','fband','f')
end
if IPso==4
  if time(1,1)==time(end,1)
      optfilename=['opti_fmin_' sname '_' num2str(time(1,1)) '.mat'];
  else
      optfilename=['opti_fmin_' sname '_' num2str(time(1,1)) '_' num2str(time(end,1)) '.mat'];
  end
  load(optfilename)   
end
if IPso==2
  if time(1,1)==time(end,1)
      optfilename=['opti_pso_' sname '_' num2str(time(1,1)) '.mat'];
  else
%      optfilename=['opti_pso_' sname '_35'  '.mat'];
      optfilename=['opti_pso_' sname '_' num2str(time(1,1)) '_' num2str(time(end,1)) '.mat'];
  end
  load(optfilename)   
end
%%
%PART IV: Harmonic analysis
% Extracting river-tide signal and capturing s(t) and f(t) as well as the time-dependent results of tidal constituents£¨Eta,Phi)
% note: the amplitude and phase of each tidal constituents is not callibrated by inference tide and nodal correction
% cof=cof';
if tsnr==1
    %error estimation with time,tidecon include timely snr£¬famj(amplitude),pha(phase) of constituent,emaj and epha are error value
    %st,ft,yout and yout-snr is reconstruct water level
 [st,ft,yout,yout_snr,percent,si,b,Eta,Phi,tidecon]=Rtide_harmonic_witherr(xin,Q,t,fu,cof,vu,f,Qc,fband,synth);
else
[st,ft,yout,Eta,Phi,percent,si,b]=Rtide_harmonic(xin,Q,t,fu,cof,vu,f,Qc,fband);
end
%--------Generate a 'prediction' using significant constituents----------
if ipre==1
    if time(1,1)==time(end,1)
         fname1=['harmcofforpre_' sname '_' num2str(time(1,1)) '.mat'];
    else
        fname1=['harmcofforpre_' sname '_' num2str(time(1,1)) '_' num2str(time(end,1)) '.mat'];
    end
    save(fname1,'cof','fu','b','Qc','fband','ju','nameu','f','vu')
%      yout= Rtide_predict(xin,Q,t,fu,cof,b,vu,f,Qc,fband);
% elseif ipre==2
%     [percent,si,percent1,si1]=Rtide_pre(xin,Q,t,fu,cof,vu,f,Qc,fband);
end

