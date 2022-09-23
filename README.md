# R_TIDE-V1.0-Matlab-Toolbox

Huayang Cai

Institute of Estuarine and Coastal Research, School of Ocean Engineering and Technology, Sun Yat-sen University, Guangzhou, 510275, China
Correspondence: Huayang Cai (caihy7@mail.sysu.edu.cn)

2022/09/23


A data-driven model to quantify the impact of river discharge on tide-river dynamics in river deltas.
Version 1.0 - September 2022

Provided by Huayang Cai

Institute of Estuarine and Coastal Research, School of Ocean Engineering and Technology, Sun Yat-sen University

Email contacts: caihy7@mail.sysu.edu.cn

How to cite:

- Huayang Cai, Bo Li, Erwan Garel, Tongtiegang Zhao, Feng Liu, and Suying Ou (2022), A data-driven model to quantify the impact of river discharge on tide-river dynamics in the Yangtze River estuary, Journal of Hydrology, submitted

In this newly proposed version, we:

1) proposed a new definition of critical river discharge Qc defined as the value that leads to an apparent shift of tidal phases by approximately 100~180° with the increase of river discharge.

2) added an error estimation model.


How to use R_TIDE

1.	Download and install R_TIDE toolbox
Users can download the latest R_TIDE toolbox from Github:
https://github.com/Huayangcai/R_TIDE-V1.0-Matlab-Toolbox.git

2.	R_TIDE Demo

2.1.	Harmonic analysis driven by river discharge

First of all, you need to load data provided by R_TIDE Toolbox (such as `Data_Yangtze_river.mat`). The demo can be executed using the main program labelled by `R_demo_Yangtze.m`.
The syntax of the main subroutine is illustrated below:

[nameu,fu,yout,st,ft,Eta,Phi,percent,si,cof]=R_tide(xin,Q,T,lat1,ray,synth,Qc,twin,sname,ipso,ipre);

Descriptions of the inputs:

`xin1`: hourly water level data used for harmonic analysis

`Q1`: hourly river discharge data used for harmonic analysis

`T1`: the corresponding time series of the input data in term of ‘datenum’

`lat1`: the latitude of the selected station

`ray`: Rayleigh criteria, the default value is 1, which indicates that the Rayleigh criteria is used to select tidal constituents. Otherwise, Rayleigh criteria is not used.

`synth`: signal noise ratio, the default value is 10. You can set it depending on your own time series.

`Qc1`: the critical discharge beyond which the tide is vanishing, the default value is the maximum of the Q variable. Generally, it can be set to be the value corresponding with a negligible tidal range last for more than 2 days.

`twin`: window spectrum, the default value is 366.

`sname`: the name of the selected station

`ipso`: the method you used for harmonic analysis. If ipso=1, it will invoke standard PSO (Particle Swarm Optimization) algorithm to optimize and save the optimized results. If ipso=2, it will directly invoke the optimized file derived from PSO last time. If ipso=3, it will invoke default Matlab FMINCON Function to optimize and save the optimized results. If ipso=4, it will directly invoke the optimized file derived from FMINCON function last time. So if it is the first time to use it, this parameter should be set as 1 or 3.

`ipre`: the method you used for prediction. If ipre=1, it will save the coefficients for prediction and then invoke the program, namely, `Rtide_predict.m`; If ipre=2, it will directly predict. Divide the time series into 2 parts, one of which is used to derive coefficient for prediction, the other is used for prediction. Users should set the length of time series in the file, namely, `Rtide_pre.m`. So if it is the first time to use it, this parameter should be set as 1.

`tsnr`: decide whether you’d like to add an error estimation. If tsnr=0, it will output the tidal properties. If tsnr=1, it will add an error estimation invoking the function 'Rtide_harmonic_witherr.m';

`tau`: the default value of travelling time of river discharge propagating to the studied tidal gauging station.


Descriptions of the outputs:

`nameu`: the name of the selected tidal constituents

`fu`: the frequency of the selected tidal constituents (/h)

`yout`: the reconstructed water level derived from R_TIDE, consisting of st and ft

`st`: the reconstructed water level derived by the residual water level model

`ft`: the reconstructed water level derived by the tidal-fluvial model

`Eta`: the time-dependent amplitude of each tidal constituent derived from R_TIDE

`Phi`: the time-dependent phase of each tidal constituent derived from R_TIDE

`precent`: the correlation coefficient between reconstructed and observed water levels

`si`: the Root Mean Square Error (RMSE) between reconstructed and observed water levels

`cof`: the regression coeffiicents adopted for each tidal bands derived from R_TIDE

`tidecon`: the tidal properties together with their errors;
