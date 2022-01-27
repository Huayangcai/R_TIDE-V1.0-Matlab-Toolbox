function [snr,tidecon]=error_snr(xres,fu,dt,errcalc,ap,am,f,vu)
% --------------Error Bar Calculations---------------------------------
%
% Error bar calcs involve two steps:
%      1) Estimate the uncertainties in the analyzed amplitude
%         for both + and - frequencies (i.e., in 'ap' and 'am').
%         A simple way of doing this is to take the variance of the
%         original time series and divide it into the amount appearing
%         in the bandwidth of the analysis (approximately 1/length).
%         A more sophisticated way is to assume "locally white"
%         noise in the vicinity of, e.g., the diurnal consistuents.
%         This takes into account slopes in the continuum spectrum.
%
%      2) Transform those uncertainties into ones suitable for ellipse
%         parameters (axis lengths, angles). This can be done 
%         analytically for large signal-to-noise ratios. However, the 
%         transformation is non-linear at lows SNR, say, less than 10
%         or so.
%

xr=fixgaps(xres); % Fill in "internal" NaNs with linearly interpolated
                  % values so we can fft things.
nreal=1;

if strmatch(errcalc(2:end),'boot'),
  fprintf('   Using nonlinear bootstrapped error estimates\n');
  
  % "noise" matrices are created with the right covariance structure
  % to add to the analyzed components to create 'nreal' REPLICATES. 
  % 

  nreal=300;             % Create noise matrices 
  [NP,NM]=noise_realizations(xr(isfinite(xr)),fu,dt,nreal,errcalc);
      
  % All replicates are then transformed (nonlinearly) into ellipse 
  % parameters.  The computed error bars are then based on the std
  % dev of the replicates.

  AP=ap(:,ones(1,nreal))+NP;        % Add to analysis (first column
  AM=am(:,ones(1,nreal))+NM;        % of NM,NP=0 so first column of
                                    % AP/M holds ap/m).
  epsp=angle(AP)*180/pi;            % Angle/magnitude form:
  epsm=angle(AM)*180/pi;
  ap=abs(AP);
  am=abs(AM);
elseif strmatch(errcalc,'linear'),
  fprintf('   Using linearized error estimates\n');
  %
  % Uncertainties in analyzed amplitudes are computed in different
  % spectral bands. Real and imaginary parts of the residual time series
  % are treated separately (no cross-covariance is assumed).
  %
  % Noise estimates are then determined from a linear analysis of errors,
  % assuming that everything is uncorrelated. This is OK for scalar time
  % series but can fail for vector time series if the noise is not 
  % isotropic.
  
  [ercx,eicx]=noise_stats(xr(isfinite(xr)),fu,dt);
  % Note - here we assume that the error in the cos and sin terms is 
  % equal, and equal to total power in the encompassing frequency bin. 
  % It seems like there should be a factor of 2 here somewhere but it 
  % only works this way! <shrug>
  [emaj,emin,einc,epha]=errell(ap+am,i*(ap-am),ercx,ercx,eicx,eicx);

  epsp=angle(ap)*180/pi;
  epsm=angle(am)*180/pi;
  ap=abs(ap);
  am=abs(am);
else
  error(['Unrecognized type of error analysis: ''' errcalc ''' specified!']);
end;

%-----Convert complex amplitudes to standard ellipse parameters--------

aap=ap./f(:,ones(1,nreal));	% Apply nodal corrections and
aam=am./f(:,ones(1,nreal));	% compute ellipse parameters.

fmaj=aap+aam;                   % major axis
fmin=aap-aam;                   % minor axis

gp=mod( vu(:,ones(1,nreal))-epsp ,360); % pos. Greenwich phase in deg.
gm=mod( vu(:,ones(1,nreal))+epsm ,360); % neg. Greenwich phase in deg.

finc= (epsp+epsm)/2;
finc(:,1)=mod( finc(:,1),180 ); % Ellipse inclination in degrees
				% (mod 180 to prevent ambiguity, i.e., 
				% we always ref. against northern 
				% semi-major axis.
	
finc=cluster(finc,180); 	% Cluster angles around the 'true' 
                                % angle to avoid 360 degree wraps.

pha=mod( gp+finc ,360); 	% Greenwich phase in degrees.

pha=cluster(pha,360);		% Cluster angles around the 'true' angle
				% to avoid 360 degree wraps.

%----------------Generate 95% CI---------------------------------------
%% For bootstrapped errors, we now compute limits of the distribution.
if strmatch(errcalc(2:end),'boot'),
     %% std dev-based estimates.
     % The 95% CI are computed from the sigmas
     % by a 1.96 fudge factor (infinite degrees of freedom).
     % emaj=1.96*std(fmaj,0,2);
     % emin=1.96*std(fmin,0,2);
     % einc=1.96*std(finc,0,2);
     % epha=1.96*std(pha ,0,2);
     %% Median-absolute-deviation (MAD) based estimates.
     % (possibly more stable?)
     % Sep 2014: I wrote 0.6375 in the next 4 lines but it should be 0.6745 
     % - thanks to Dan Codiga  for pointing this out (and Jacob  Wenegrat for 
     %   finding it again...)
      emaj=median(abs(fmaj-median(fmaj,2)*ones(1,nreal)),2)/.6745*1.96;
      emin=median(abs(fmin-median(fmin,2)*ones(1,nreal)),2)/.6745*1.96;
      einc=median(abs(finc-median(finc,2)*ones(1,nreal)),2)/.6745*1.96;
      epha=median(abs( pha-median( pha,2)*ones(1,nreal)),2)/.6745*1.96;
else
   % In the linear analysis, the 95% CI are computed from the sigmas
   % by this fudge factor (infinite degrees of freedom).
   emaj=1.96*emaj;
   emin=1.96*emin;
   einc=1.96*einc;
   epha=1.96*epha;
end;
				  
    tidecon=[fmaj(:,1),emaj,pha(:,1),epha];

% % Sort results by frequency (needed if anything has been inferred since 
% % these are stuck at the end of the list by code above).
%  [fu,I]=sort(fu);
%  nameu=nameu(I,:);
%  tidecon=tidecon(I,:);

snr=(tidecon(:,1)./tidecon(:,2)).^2;  % signal to noise ratio
