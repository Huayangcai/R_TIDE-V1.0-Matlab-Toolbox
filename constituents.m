%----------------------------------------------------------------------
function [nameu,fu,ju,namei,fi,jinf,jref]=constituents(minres,constit,...
                                     shallow,infname,infref,centraltime);
% [name,freq,kmpr]=constituents(minres,infname) loads tidal constituent
% table (containing 146 constituents), then picks out only the '
% resolvable' frequencies (i.e. those that are MINRES apart), base on 
% the comparisons in the third column of constituents.dat. Only 
% frequencies in the 'standard' set of 69 frequencies are actually used.
% Also return the indices of constituents to be inferred.

% If we have the mat-file, read it in, otherwise create it and read
% it in!

% R Pawlowicz 9/1/01 
% Version 1.0
%
%    19/1/02 - typo fixed (thanks to  Zhigang Xu)

% Compute frequencies from astronomical considerations.

if minres>1/(18.6*365.25*24),                       % Choose only resolveable pairs for short
  [const,sat,cshallow]=t_getconsts(centraltime);    % Time series   
  ju=find(const.df>=minres);
else                                                % Choose them all if > 18.6 years.
  [const,sat,cshallow]=t_get18consts(centraltime);
  ju=[2:length(const.freq)]';  % Skip Z0
  for ff=1:2,                  % loop twice to make sure of neightbouring pairs
    jck=find(diff(const.freq(ju))<minres);
    if (length(jck)>0)
       jrm=jck;
       jrm=jrm+(abs(const.doodsonamp(ju(jck+1)))<abs(const.doodsonamp(ju(jck))));
       disp(['  Warning! Following constituent pairs violate Rayleigh criterion']);
       for ick=1:length(jck);
	 disp(['     ',const.name(ju(jck(ick)),:),' vs ',const.name(ju(jck(ick)+1),:) ' - not using ',const.name(ju(jrm(ick)),:)]);
       end;
       ju(jrm)=[];
    end
  end;
end;
  
if ~isempty(constit),     % Selected if constituents are specified in input.
  ju=[];
  for k=1:size(constit,1),
   j1=strmatch(constit(k,:),const.name);
   if isempty(j1),
     disp(['Can''t recognize name ' constit(k,:) ' for forced search']);
   elseif j1==1,
     disp(['*************************************************************************']);
     disp(['Z0 specification ignored - for non-tidal offsets see ''secular'' property']);
     disp(['*************************************************************************']);
   else  
     ju=[ju;j1];
   end;
  end;
  [dum,II]=sort(const.freq(ju)); % sort in ascending order of frequency.
  ju=ju(II);
end;


disp(['   number of standard constituents used: ',int2str(length(ju))])

if ~isempty(shallow),          % Add explictly selected shallow water constituents.
 for k=1:size(shallow,1),
   j1=strmatch(shallow(k,:),const.name);
   if isempty(j1),
     disp(['Can''t recognize name ' shallow(k,:) ' for forced search']);
   else
     if isnan(const.ishallow(j1)),
       disp([shallow(k,:) ' Not a shallow-water constituent']);
     end;
     disp(['   Forced fit to ' shallow(k,:)]);
     ju=[ju;j1];
   end;
 end;
 
end;
      
nameu=const.name(ju,:);
fu=const.freq(ju);


% Check if neighboring chosen constituents violate Rayleigh criteria.
jck=find(diff(fu)<minres);
if (length(jck)>0)
   disp(['  Warning! Following constituent pairs violate Rayleigh criterion']);
   for ick=1:length(jck);
   disp(['     ',nameu(jck(ick),:),'  ',nameu(jck(ick)+1,:)]);
   end;
end

% For inference, add in list of components to be inferred.

fi=[];namei=[];jinf=[];jref=[];
if ~isempty(infname),
  fi=zeros(size(infname,1),1);
  namei=zeros(size(infname,1),4);
  jinf=zeros(size(infname,1),1)+NaN;
  jref=zeros(size(infname,1),1)+NaN;

  for k=1:size(infname,1),
   j1=strmatch(infname(k,:),const.name);
   if isempty(j1),
     disp(['Can''t recognize name' infname(k,:) ' for inference']);
   else
    jinf(k)=j1;
    fi(k)=const.freq(j1);
    namei(k,:)=const.name(j1,:);
    j1=strmatch(infref(k,:),nameu);
    if isempty(j1),
      disp(['Can''t recognize name ' infref(k,:) ' for as a reference for inference']);
    else
      jref(k)=j1;
      fprintf(['   Inference of ' namei(k,:) ' using ' nameu(j1,:) '\n']);
    end;
   end;
  end;    
  jinf(isnan(jref))=NaN;
end;
