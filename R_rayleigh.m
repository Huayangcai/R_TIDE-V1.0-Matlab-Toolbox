function minres = R_rayleigh(t,Q,xtune,fu,win,crit,minres)

%The function ns_rayleigh calculates a new rayleigh criterion based on 
%the spectral content of each term of the basis functions and returns the
%frequencies that meet this criterion.
%INPUT: t:          time vector
%       Q:          flow vector/matrix
%       R:          tidal range vector/matrix
%       P:          atmospheric vector/matrix
%       xtune:      exponents in the basis functions
%       Qc:         discharge cutoff values (optional)
%       fu:         initial frequencies (model-specific)
%       win:        window for the PSD analysis (should be one year)
%       crit:       fraction of the total power on which is based the criterion
%                   (1-crit) => cumulative integral from 0 to w
%       minres:     standard rayleigh criterion
%       df:         frequency difference based on the specified decision tree
%       addfu:      added frequencies
%The other input arguments are passed to ns_get_constituents function.

%Terms of the basis functions
%Here we use the stage model with the exponents of the tidal-fluvial model 
%only to get the intensity of the modulation for each non-stationary term

% w_lim = 0.001; %limit of the first cusp

[mQ nQ] = size(Q);

timestep = round((t(end)-t(1))/length(t)*24*60)/60; %in hours
tc = []; 
tc_aux = [];

          tc = [tc, Q(:,1).^xtune];
          tc = [tc, cos((2*pi*t)*fu'), sin((2*pi*t)*fu')];
            
            %discharge term(s)
            for i=1:nQ
                tc = [tc, (Q(:,i).^xtune)*ones(1,length(fu)).*cos((2*pi*t)*fu'),...
                          (Q(:,i).^xtune)*ones(1,length(fu)).*sin((2*pi*t)*fu')];
            end
  
%[ms,ns]=size(tc)
%w_cr=zeros(nQ,1);
for i=1:nQ %loop on the total number of terms
    term = tc(:,i);  term = term(~isnan(term));
    term = term - mean(term);
    win  = min(win,length(term));
    [Pf w] = pwelch(term,win); %in units of power per radians per sample
    w = w/2/pi/timestep;
%     ilim = find(w<=w_lim,1,'last');
%     w = w(1:ilim);
%     Pf = Pf(1:ilim);
    %in units of cycles per sample
%     [Pf w] = psd(term,win,1); %one can use PWELCH instead but results are sligthly different
    %Cumulative integrals
    cumP = cumtrapz(w,Pf);
    %Normalizing
    cumP = cumP/cumP(end);
    %Frequency corresponding to (1-crit)% of the surface
    k = find(cumP>1-crit,1,'first');
    if ~isempty(k), w_cr(i) = w(k-1) + (w(k)-w(k-1))/(cumP(k)-cumP(k-1))*(1-crit - cumP(k-1)); end
end

%Rayleigh criterion
minres_cr = max(w_cr);
minres = max([minres; minres_cr]);

% minres = opt.ray/(max(data.t)-min(data.t));




