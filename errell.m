%----------------------------------------------------------------------
function [emaj,emin,einc,epha]=errell(cxi,sxi,ercx,ersx,ercy,ersy)
% [emaj,emin,einc,epha]=errell(cx,sx,cy,sy,ercx,ersx,ercy,ersy) computes
% the uncertainities in the ellipse parameters based on the 
% uncertainities in the least square fit cos,sin coefficients.
%
%  INPUT:  cx,sx=cos,sin coefficients for x 
%          cy,sy=cos,sin coefficients for y
%          ercx,ersx=errors in x cos,sin coefficients
%          ercy,ersy=errors in y cos,sin coefficients
%          
%  OUTPUT: emaj=major axis error
%          emin=minor axis error
%          einc=inclination error (deg)
%          epha=pha error (deg)

% based on linear error propagation, with errors in the coefficients 
% cx,sx,cy,sy uncorrelated. 

% B. Beardsley  1/15/99; 1/20/99
% Version 1.0

r2d=180./pi;
cx=real(cxi(:));sx=real(sxi(:));cy=imag(cxi(:));sy=imag(sxi(:));
ercx=ercx(:);ersx=ersx(:);ercy=ercy(:);ersy=ersy(:);

rp=.5.*sqrt((cx+sy).^2+(cy-sx).^2);
rm=.5.*sqrt((cx-sy).^2+(cy+sx).^2);
ercx2=ercx.^2;ersx2=ersx.^2;
ercy2=ercy.^2;ersy2=ersy.^2;

% major axis error
ex=(cx+sy)./rp;
fx=(cx-sy)./rm;
gx=(sx-cy)./rp;
hx=(sx+cy)./rm;
dcx2=(.25.*(ex+fx)).^2;
dsx2=(.25.*(gx+hx)).^2;
dcy2=(.25.*(hx-gx)).^2;
dsy2=(.25.*(ex-fx)).^2;
emaj=sqrt(dcx2.*ercx2+dsx2.*ersx2+dcy2.*ercy2+dsy2.*ersy2);

% minor axis error
dcx2=(.25.*(ex-fx)).^2;
dsx2=(.25.*(gx-hx)).^2;
dcy2=(.25.*(hx+gx)).^2;
dsy2=(.25.*(ex+fx)).^2;
emin=sqrt(dcx2.*ercx2+dsx2.*ersx2+dcy2.*ercy2+dsy2.*ersy2);

% inclination error
rn=2.*(cx.*cy+sx.*sy);
rd=cx.^2+sx.^2-(cy.^2+sy.^2);
den=rn.^2+rd.^2;
dcx2=((rd.*cy-rn.*cx)./den).^2;
dsx2=((rd.*sy-rn.*sx)./den).^2;
dcy2=((rd.*cx+rn.*cy)./den).^2;
dsy2=((rd.*sx+rn.*sy)./den).^2;
einc=r2d.*sqrt(dcx2.*ercx2+dsx2.*ersx2+dcy2.*ercy2+dsy2.*ersy2);

% phase error
rn=2.*(cx.*sx+cy.*sy);
rd=cx.^2-sx.^2+cy.^2-sy.^2;
den=rn.^2+rd.^2;
dcx2=((rd.*sx-rn.*cx)./den).^2;
dsx2=((rd.*cx+rn.*sx)./den).^2;
dcy2=((rd.*sy-rn.*cy)./den).^2;
dsy2=((rd.*cy+rn.*sy)./den).^2;
epha=r2d.*sqrt(dcx2.*ercx2+dsx2.*ersx2+dcy2.*ercy2+dsy2.*ersy2);