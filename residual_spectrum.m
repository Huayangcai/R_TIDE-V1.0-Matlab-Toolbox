%----------------------------------------------------------------------
function [fband,Pxrave,Pxiave,Pxcave]=residual_spectrum(xres,fu,dt)
% RESIDUAL_SPECTRUM: Computes statistics from an input spectrum over
% a number of bands, returning the band limits and the estimates for
% power spectra for real and imaginary parts and the cross-spectrum.          
%
% Mean values of the noise spectrum are computed for the following 
% 8 frequency bands defined by their center frequency and band width:
% M0 +.1 cpd; M1 +-.2 cpd; M2 +-.2 cpd; M3 +-.2 cpd; M4 +-.2 cpd; 
% M5 +-.2 cpd; M6 +-.21 cpd; M7 (.26-.29 cpd); and M8 (.30-.50 cpd). 

% S. Lentz  10/28/99
% R. Pawlowicz 11/1/00
% Version 1.0

% Define frequency bands for spectral averaging.
fband =[.00010 .00417;
        .03192 .04859;
        .07218 .08884;
        .11243 .12910;
        .15269 .16936;
        .19295 .20961;
        .23320 .25100;
        .26000 .29000;
        .30000 .50000];

% If we have a sampling interval> 1 hour, we might have to get
% rid of some bins.
%fband(fband(:,1)>1/(2*dt),:)=[];

nfband=size(fband,1);
nx=length(xres);

% Spectral estimate (takes real time series only).


% Matlab has changed their spectral estimator functions
% To match the old code, I have to divide by 2*dt. This is because
% 
%  PSD*dt  is two-sided spectrum in units of power per hertz.
%
%  PWELCH is the one-sided spectrum in power per hertz
%
%  So PWELCH/2 = PSD*dt


%[Pxr,fx]=psd(real(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
%[Pxi,fx]=psd(imag(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
%[Pxc,fx]=csd(real(xres),imag(xres),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.


[Pxr,fx]=pwelch(real(xres),hanning(nx),ceil(nx/2),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme. If you have an error here you are probably missing this toolbox
Pxr=Pxr/2/dt;
[Pxi,fx]=pwelch(imag(xres),hanning(nx),ceil(nx/2),nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
Pxi=Pxi/2/dt;
[Pxc,fx]=cpsd(real(xres),imag(xres),[],[],nx,1/dt); % Call to SIGNAL PROCESSING TOOLBOX - see note in t_readme.
Pxc=Pxc/2/dt;

df=fx(3)-fx(2);
Pxr(round(fu./df)+1)=NaN ; % Sets Px=NaN in bins close to analyzed frequencies
Pxi(round(fu./df)+1)=NaN ; % (to prevent leakage problems?).
Pxc(round(fu./df)+1)=NaN ; 

Pxrave=zeros(nfband,1);
Pxiave=zeros(nfband,1);
Pxcave=zeros(nfband,1);
% Loop downwards in frequency through bands (cures short time series
% problem with no data in lowest band).
%
% Divide by nx to get power per frequency bin, and multiply by 2
% to account for positive and negative frequencies.
%
for k=nfband:-1:1,
   jband=find(fx>=fband(k,1) & fx<=fband(k,2) & isfinite(Pxr));
   if any(jband),
     Pxrave(k)=mean(Pxr(jband))*2/nx;
     Pxiave(k)=mean(Pxi(jband))*2/nx;
     Pxcave(k)=mean(Pxc(jband))*2/nx;
   elseif k<nfband,
     Pxrave(k)=Pxrave(k+1);   % Low frequency bin might not have any points...
     Pxiave(k)=Pxiave(k+1);   
     Pxcave(k)=Pxcave(k+1);   
   end;
end
 

