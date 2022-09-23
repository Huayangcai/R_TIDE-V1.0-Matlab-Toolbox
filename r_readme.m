% r_readme.m
%
%                River-Tidal Analysis Toolbox R_TIDE
%                              by        
%               Huayang Cai (Sun Yat-sen University) 
%
%                 Version 2.0    September 21, 2022
%
%Please citing:
%1. Huayang Cai, Bo Li, Erwan Garel, Tongtiegang Zhao, Feng Liu, and Suying Ou (2022), A data-driven model to quantify the impact of river discharge on tide-river dynamics in the Yangtze River estuary, Journal of Hydrology, submitted
%
% The toolbox presently contains the following mfiles:
%
% ---FOR ANALYSIS
% r_tide.m       - computes the tidal analysis of the real tidal 
%                  time series using River-tidal Data-driven Model (R_TIDE).
% cluster.m   -Cluster angles in rows around the angles in the first column.
%
% t_vuf.m        - computes nodal corrections (from T_TIDE package). 
% 
% t_astron.m     - computes astronomical arguments (from T_TIDE package). 
%
% t_getconsts.m     - Gets constituent data structures (from T_TIDE package).
%
% t_get18consts.m     - the same as t_getconsts.m.
%
% constituents.m     - Compute frequencies from astronomical considerations.
%
% errell.m     - compute the uncertainities in the ellipse parameters based on the uncertainities in the least square fit cos,sin coefficients.
%
% error_snr.m     - Error Bar Calculations.
%
% fixgaps.m     - Linearly interpolate gaps in a time series
%
% fmin_Q_tide.m     - derive Qc via FMINCON method
%
% noise_realizations.m     - Generates matrices of noise (with correct cross-correlation structure) for bootstrap analysis.
%
% noise_stats.m     - Computes statistics of residual energy for all constituents
%
% Qztidefun.m     - Computes the minimum RMSE for the model
%
% R_Qtidal.m     - Computes the critical river discharge
%
% R_rayleigh.m     - Rayleigh Criteria
%
% residual_spectrum.m     - Computes statistics from an input spectrum over a number of bands, returning the band limits and the estimates for power spectra for real and imaginary parts and the cross-spectrum. 
%
% Rtide_discrete.m     - discrete diurnal and semi-diurnal signal
%
% Rtide_harmonic.m     - non-stationary harmonic analysis by the optimized cof and the driving Q.
%
% Rtide_harmonic_witherr.m     - non-stationary harmonic analysis by the optimized cof and the driving Q and its error estimation.
%
% Rtide_pre.m     - derive the optimized coefficients for prediction.
%
% Rtide_predict.m     - predict water level by the cofficent b of sine and cosine function, the optimized cof and the driving Q  
%
% Rtide_pso.m     - Particle Swarm Optimation algorithm 
%
% Rtide_sixcof.m     - least square method for R2 or RMSE
%
% ---FOR DOCUMENTATION
%
% r_readme.m     - this file.
%
% ---FOR DEMONSTRATION
%
% R_demo_Yangtze.m   - a short example using the Yangtze River Estuary water level data.
%
% Various data files are also included:
%
% t_constituents.mat - constituent data structures.
%
% constituents.mat - constituent data structures.
%
% t_equilib.mat - Equilibrium Tidal Potential (Values taken from Appendix 1 of Godin, "The Analysis of Tides", Univ. Toronto Press, 1972.)
%
% tide3.dat      - standard constituent data file from the IOS analysis package.
%
%
%
% Questions or comments to:
% Huayang Cai (Sun Yat-sen University) 
% Email:caihy7@mail.sysu.edu.cn

help r_readme