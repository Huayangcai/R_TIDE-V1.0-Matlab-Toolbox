function [NP,NM]=noise_realizations(xres,fu,dt,nreal,errcalc);
% NOISE_REALIZATIONS: Generates matrices of noise (with correct
% cross-correlation structure) for bootstrap analysis.
%

% R. Pawlowicz 11/10/00
% Version 1.0

if strmatch(errcalc,'cboot'),
  [fband,Pxrave,Pxiave,Pxcave]=residual_spectrum(xres,fu,dt);
  
  Pxcave=zeros(size(Pxcave));  %% For comparison with other technique!
  %fprintf('**** Assuming no covariance between u and v errors!*******\n');

elseif strmatch(errcalc,'wboot'),
  fband=[0 .5];
  nx=length(xres);
  A=cov(real(xres),imag(xres))/nx;
  Pxrave=A(1,1);Pxiave=A(2,2);Pxcave=A(1,2);
else
  error(['Unrecognized type of bootstap analysis specified: ''' errcalc '''']);
end;
  
nfband=size(fband,1);

Mat=zeros(4,4,nfband);
for k=1:nfband,

  % The B matrix represents the covariance matrix for the vector
  % [Re{ap} Im{ap} Re{am} Im{am}]' where Re{} and Im{} are real and
  % imaginary parts, and ap/m represent the complex constituent 
  % amplitudes for positive and negative frequencies when the input
  % is bivariate white noise. For a flat residual spectrum this works 
  % fine.
 
  % This is adapted here for "locally white" conditions, but I'm still
  % not sure how to handle a complex sxy, so this is set to zero
  % right now.
  
  p=(Pxrave(k)+Pxiave(k))/2;
  d=(Pxrave(k)-Pxiave(k))/2;
  sxy=Pxcave(k);
  
  B=[p    0   d   sxy;
     0    p  sxy  -d;
     d   sxy  p    0
     sxy -d   0    p];

  % Compute the transformation matrix that takes uncorrelated white 
  % noise and makes noise with the same statistical structure as the 
  % Fourier transformed noise.
  [V,D]=eig(B);
  Mat(:,:,k)=V*diag(sqrt(diag(D)));
end;

% Generate realizations for the different analyzed constituents.

N=zeros(4,nreal);
NM=zeros(length(fu),nreal);
NP=NM;
for k=1:length(fu);
  l=find(fu(k)>fband(:,1) & fu(k)<fband(:,2));
  N=[zeros(4,1),Mat(:,:,l)*randn(4,nreal-1)];
  NP(k,:)=N(1,:)+i*N(2,:);
  NM(k,:)=N(3,:)+i*N(4,:);
end;
