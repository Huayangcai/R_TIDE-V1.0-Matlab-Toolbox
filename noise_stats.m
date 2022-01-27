

%----------------------------------------------------------------------
function [ercx,eicx]=noise_stats(xres,fu,dt);
% NOISE_STATS: Computes statistics of residual energy for all 
% constituents (ignoring any cross-correlations between real and
% imaginary parts).

% S. Lentz  10/28/99
% R. Pawlowicz 11/1/00
% Version 1.0

[fband,Pxrave,Pxiave,Pxcave]=residual_spectrum(xres,fu,dt);
nfband=size(fband,1);
mu=length(fu);

% Get the statistics for each component.
ercx=zeros(mu,1);
eicx=zeros(mu,1);
for k1=1:nfband;
   k=find(fu>=fband(k1,1) & fu<=fband(k1,2));
   ercx(k)=sqrt(Pxrave(k1));
   eicx(k)=sqrt(Pxiave(k1));
end