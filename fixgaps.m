%----------------------------------------------------------------------
function y=fixgaps(x);
% FIXGAPS: Linearly interpolates gaps in a time series
% YOUT=FIXGAPS(YIN) linearly interpolates over NaN in the input time 
% series (may be complex), but ignores trailing and leading NaNs.

% R. Pawlowicz 11/6/99
% Version 1.0

y=x;

bd=isnan(x);
gd=find(~bd);

bd([1:(min(gd)-1) (max(gd)+1):end])=0;


y(bd)=interp1(gd,x(gd),find(bd)); 


%----------------------------------------------------------------------


