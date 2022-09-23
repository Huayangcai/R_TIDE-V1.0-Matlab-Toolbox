% R_DEMO - demonstration of capabilities used in Yangtze River Estuary,

% The parameters input:
% xin1 -- water level data used for harmonic analysis
% Q1 -- discharge data used for harmonic analysis
% T1 -- the corresponding time series of the input data in term of 'datenum'
% lat1 -- the latitude of the selected station
% ray -- Rayleigh criteria, the default value is 1, which indicates that the Rayleigh criteria is used to select tidal constituents. Otherwise, Rayleigh criteria is not used.
% synth -- signal noise ratio, the default value is 10. You can set it depending on your own time series
% Qc1 -- the critical discharge beyond which the tide is vanishing, the default value is the maximum of the Q variable.  Generally, it can be set to be the value corresponding with a negligible tidal range last for more than 2 days
% twin -- window spectrum, the default value is 366.
% sname -- the name of the selected station
% ipso -- the method you used for harmonic analysis. If ipso=1, it will invoke PSO to optimize and save the optimized results. If ipso=2, it will directly invoke the optimized file derived last time. So if it is the first time to use it, this parameter should be set as 1.
% ipre -- the method you used for prediction. If ipre=1, it will save the coefficients for prediction and then invoke the program, namely, Rtide_predict.m; If ipre=2, it will directly predict. Divide the data series into 2 parts, one of which is used to derive coefficient for prediction, another one is used to predict. Users should set the length of data series in the file, namely, Rtide_pre.m. So if it is the first time to use it, this parameter should be set as 1.
% tsnr -- error estimation. if tsnr=0, derive tidal properties; if tsnr=1, add an error estimation.
% tau -- the default value of travelling time;

% The parameters output:
% nameu -- the name of the selected tidal constituents
% fu -- the frequency of the selected tidal constituents (/h)
% yout -- the reconstructed water level derived from R_TIDE, consisting of st and ft
% st -- the reconstructed water level derived by the residual water level model
% ft -- the reconstructed water level derived by the tidal-fluvial model
% Eta -- the time-dependent amplitude of each tidal constituent derived from R_TIDE
% Phi -- the time-dependent phase of each tidal constituent derived from R_TIDE
% precent -- the correlation coefficient between reconstructed and observed water levels
% si -- the Root Mean Square Error (RMSE) between reconstructed and observed water levels
% cof -- the regression coefficients adopted for each tidal constituent bands derived from R_TIDE
% tidecon -- the tidal properties together with their errors;

% Reference:
% Contract: Huayang Cai (E-mail: caihy7@mail.sysu.edu.cn); Suying Ou (E-mail: ousuying@mail.sysu.edu.cn)
%%
clc,clear
close all

load Data_Yangtze_river.mat   % Load the example.
mm=find(ZQ(:,1)==datenum(2009,1,1,0,0,0));
ii1=1:mm-1; % select the time series for calibration (2002-2008)
T1=ZQ(ii1,1); % time series from 2002 to 2008
Q1=ZQ(ii1,end); % river discharge from 2002 to 2008
Z1=ZQ(ii1,2:end-1); % water levels from 2002 to 2008
ii2=mm:length(ZQ);%select time series for prediction (2009-2012);
T2=ZQ(ii2,1); % time series from 2002 to 2008
Q2=ZQ(ii2,end); % river discharge from 2002 to 2008
Z2=ZQ(ii2,2:end-1); % water levels from 2002 to 2008
Qc1=max(Q1); % the default value of critical river discharge
ray=1; % Rayleigh criteria
synth=1; % signal-to-noise ratio
twin=366; % days for spectrum window
tsnr=0;
% if tsnr=0, derive tidal properties tidecon;
% if tsnr=1, add the error estimation together with the output tidecon;
ipso=2; 
% if ipso=1, using PSO to optimize and save the optimized results.
% if ipso=2, directly using the optimized file derived from PSO last time.
% if ipso=3, using Fmincon function to optimize and save the optimized results.
% if ipso=4, directly using the optimized file derived from Fmincon function last time.
ipre=1;
% if ipre=1, save the coefficients for prediction (ipre=2). 
% if ipre=2, directly predict. Divide the time series into 2 parts, one of which is used to derive coefficient, the one is used for prediction.
% Users should set the length of time series in Rtide_pre.m.
ic=1;
% if ic=0, ignore the influence of large river discharge on tidal properties;
% if ic=1, need a given critical Qc, the default value is the maximum value,
% then calibrate them according to the results derived from the first test
Qc1(1:6,1)=66000;
Qc1(1:6,2)=66000;

if ic==1
    Qc1(1:6,1)=[66000;66000;66000;66000;66000;66000];
else
    Qc1(1:6,2)=max(Q1); 
end   
tau=[15;13;9;7;5;4]; % the default value of tau, repersenting the travel time of the water particle

% for i=1:size(Z1,2)
%     xin1=Z1(:,i); % read the water level of selected time series in each tidal gauging station
%     dH=0.05;% The cofficient of critical tidal range,  the corresponding river discharge is the critical discharge Qc2. If the river discharge exceeds this value, there is no tidal signal in this station.    
%     Qc0=R_Qtidal(xin1,Q1,T1,dH); % calculate the critical river discharge Qc2,but for many staion,it have no 
%     Qc1(i,1)=Qc0;%mean no second critical ;
% end
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
    [nameu,fu,yout,st,ft,Eta,Phi,percent,si,cof,tidecon]=R_tide(xin1,Q1,T1,lat1,ray,synth,Qc1(i,:),twin,sname,ipso,ipre,tsnr,tau(i));% run R_tide program
    per(i)=percent;rmse(i)=si;
    if (ipso==1 || ipso==2)
        fname1=['R_TIDE_calibration_for_PSO_' sname '.mat'];
    end
    if (ipso==3 || ipso==4)
        fname1=['R_TIDE_calibration_for_fmin_' sname '.mat'];
    end
    save(fname1,'T1','xin1','Q1','st','ft','percent','Eta','Phi','cof','si','yout','Qc1','nameu','tidecon'); % save the results   
 
 %discrete case
 %[y1,y2,yok,yms]= Rtide_discrete(xin1,Q1,T1,sname)

    %% validation
    % In this part, we use the coefficients derived from calibration period
    % as well as the river discharge records in the rest time series to
    % reconstruct the water level series, which reveal the influence of the
    % alteration in river discharge, ignoring the geometric change
 % ***************************************************************************************** 
 %prediction 
 if ipre==2
     [z1,y2,percent1,si1,st,ft]=Rtide_pre(Z2(:,i),Q2,T2,sname);
    if (ipso==1 || ipso==2)
        fname2=['R_TIDE_outputforper_for_PSO_' sname '.mat'];
    end
    if (ipso==3 || ipso==4)
        fname2=['R_TIDE_outputforper_for_fmin_' sname '.mat'];
    end
    save (fname2,'z1','y2','cof','percent1','st','ft','Q2','si1')
 end
end
