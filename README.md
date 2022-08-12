# R_TIDE-V1.0-Matlab-Toolbox

Huayang Cai

Institute of Estuarine and Coastal Research, School of Ocean Engineering and Technology, Sun Yat-sen University, Guangzhou, 510275, China
Correspondence: Huayang Cai (caihy7@mail.sysu.edu.cn)

2022/01/28


A data-driven model to quantify the impact of river discharge on tide-river dynamics in river deltas.
Version 1.0 - January 2022

Provided by Huayang Cai

Institute of Estuarine and Coastal Research, School of Ocean Engineering and Technology, Sun Yat-sen University

Email contacts: caihy7@mail.sysu.edu.cn

How to cite:

- Huayang Cai, Bo Li, Erwan Garel, Tongtiegang Zhao, Feng Liu, and Suying Ou (2022), A data-driven model to quantify the impact of river discharge on tide-river dynamics in the Yangtze River estuary, Journal of Hydrology, submitted

How to use R_TIDE

1.	Download and install R_TIDE toolbox
Users can download the latest R_TIDE toolbox from Github:
https://github.com/Huayangcai/R_TIDE-V1.0-Matlab-Toolbox.git

2.	R_TIDE Demo

2.1.	Harmonic analysis driven by river discharge

First of all, you need to load data provided by R_TIDE Toolbox (such as `Data_Yangtze_river.mat`). The demo can be executed using the main program labelled by `R_demo_Yangtze.m`.

The data file `Data_Yangtze_river.mat` contains 2 variables, including `stname` and `ZQ`. 

`stname` denotes the name of tidal gauging stations, including 6 columns (e.g., TSG, JY, ZJ, NJ, MAS, WH, respectively). 
For instance:

'TSG'

'JY'

'ZJ'

'NJ'

'MAS'

'WH'

`ZQ` denotes hourly data used for harmonic analysis. The data in the 1st column denote the time series of the input data in term of ‘datenum’. The data between the 2nd and the 7th column denote the water level series observed in the tidal stations mentioned above. For instance, there are 6 columns of water levels in this variable, the data in the 2nd column represent the water levels in TSG and the data in the 7th column represent the water levels in WH. The data in the 8th column denote hourly river discharge data used for harmonic analysis.

For instance:

731217.166666667	0.7818	-0.0565	1.4980	2.1796	1.9357	2.2544	12731.5430

731217.208333333	1.2703	0.6749	0.3553	1.3052	1.8316	2.1548	12770.1389

731217.250000000	1.4570	1.2471	0.1852	0.8740	1.7609	2.0815	12806.2174

731217.291666667	1.3719	1.5003	0.6139	0.7844	1.7172	2.0306	12839.8438

731217.333333333	1.1476	1.4408	1.2681	0.9349	1.6940	1.9980	12871.0829

……




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

