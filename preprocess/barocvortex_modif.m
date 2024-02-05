%
% Compute zeta,ubar,vbar,u,v,t for the vortex
%
%
% Pierrick Penven, IRD, 2005
% 2014: addition gr

%
%
% P1 : pressure at z=0
%
% P1=P0*exp(-r2/lambda2)
%
% P0 such as u(r0)=umax
%
% => P0=rho0 f0 umax lambda sqrt(e/2)
%
P0=rho0*f0*umax*Rmax*exp(1/alpha);
rX=sqrt(X.^2+Y.^2);
P1=Pa+P0*exp(-(rX.^alpha)/alpha/Rmax^alpha); 
%P1=Pa+P0*exp(-(X.^alpha+Y.^alpha)/alpha/Rmax^alpha);
%
% Calcul de rho at z=0
%
%sigma0=-P0*(1-exp(-H))/(g*(H-1+exp(-H)));
sigma0=-P0/(g*H)
rho1=rho0+sigma0*exp(-(rX.^alpha)/alpha/Rmax^alpha);
%
% Surface elevation 
%
zeta=(P1-Pa)./(g.*rho1);
%
% Vertical grid 
%
zw=zlevs(h0,zeta,theta_s,theta_b,hc,N,'w',vtransform);
zr=zlevs(h0,zeta,theta_s,theta_b,hc,N,'r',vtransform);
%
% Density
%
xr=tridim(X,N);
yr=tridim(Y,N);
r=sqrt(xr.^2+yr.^2);
%
%rho=rho0.*(1-N2.*zr./g);
R0=26;
T0=25;
TCOEF=0.28;
%%%%%%%%%% case if temp back defined :
%TempBack=T1-(T1-Tsurf)*exp(zr./ZT)+St*zr;
%rhoBack=1000+R0-TCOEF*(TempBack-T0);
%%%%%%%%%%% case if rho back defined :
rhoBack=1000+R1 - (R1-Rho_0)*exp(zr./ZT) - Sl*zr;

%rhodyn=-P0.*(1-exp(-zr-H)).*exp(-(xr.^alpha+yr.^alpha)/alpha/Rmax^alpha)./(g*H)
rhodyn=P0.*zr./H.*exp(-zr.^2/2/H^2).*exp(-r.^alpha/alpha/Rmax^alpha)./(g*H);
%rho(zr>-H)=rho(zr>-H)+rhodyn(zr>-H);
%rho=rhoBack;
%rho(zr>-H)=rhoBack(zr>-H)+rhodyn(zr>-H);
rho=rhoBack+rhodyn;
%
% Temperature
%
% rho=rho0+R0-TCOEF*T
%
t=T0+(-rho+1000+R0)/TCOEF;
%
% U and V
%
%a=2*P0/(f0*rho0*lambda^2);
%zu=0.5*(zr(:,:,1:end-1)+zr(:,:,2:end));
%xu=0.5*(xr(:,:,1:end-1)+xr(:,:,2:end));
%yu=0.5*(yr(:,:,1:end-1)+yr(:,:,2:end));
%F=(H-1+zu+exp(-zu-H))./(H-1+exp(-H));
%F(zu<-H)=0;
%u=a.*F.*yu.*exp(-(xu.^2+yu.^2)/lambda^2);
%zv=0.5*(zr(:,1:end-1,:)+zr(:,2:end,:));
%xv=0.5*(xr(:,1:end-1,:)+xr(:,2:end,:));
%yv=0.5*(yr(:,1:end-1,:)+yr(:,2:end,:));
%F=(H-1+zv+exp(-zv-H))./(H-1+exp(-H));
%F(zv<-H)=0;
%v=-a.*F.*xv.*exp(-(xv.^2+yv.^2)/lambda^2);
%
% Compute geostrophic velocities Vg
%
%F=(H-1+zr+exp(-zr-H))./(H-1+exp(-H));
%F(zr<-H)=0;
F=exp(-zr.^2/2/H^2);
Vg=-(r.^(alpha-1).*P0.*F./(rho0.*f0.*Rmax^alpha)).*exp(-r.^alpha/alpha/Rmax^alpha);
%
% Compute gradient wind velocities Vgr : Vgr/R^2 + f Vgr = f Vg 
%
a=1+4.*Vg./(f0.*r);
[indx]=find(a<0);
if length(indx)>0
  disp([num2str(length(indx)),...
  ' points with no gradient wind solution (use geostrophy) '])
end
a(a<0)=1;
Vgr=2.*Vg./(1+sqrt(a));

%%% Adding radial perturbance

Vr=Perturb.*Vgr;%.*sin(2*pi*CycleVert/H0.*zr);
%
%Other modes
for k=[2:ModeMax]
    RAND=rand(1);
    Vr=Vr+Perturb.*Vgr.*cos(k.*theta+RAND*2*pi);
    %Vr=Vr+Perturb.*Vtc.*cos(k.*theta+RAND*2*pi).*cos(2*pi*CycleVert/H0.*zr+RAND*2*pi);
end

%Vgr=Vg;
%
% Project on the grid
%
ur=-Vgr.*(yr)./r + Vr.*(xr)./r;
vr=Vgr.*(xr)./r + Vr.*(yr)./r;
u=0.5*(ur(:,:,1:end-1)+ur(:,:,2:end));
v=0.5*(vr(:,1:end-1,:)+vr(:,2:end,:));
%
% Barotropic speeds
%
dz=zw(2:end,:,:)-zw(1:end-1,:,:);
dzu=0.5*(dz(:,:,1:end-1)+dz(:,:,2:end));
dzv=0.5*(dz(:,1:end-1,:)+dz(:,2:end,:));
hu=squeeze(sum(dzu.*u));
hv=squeeze(sum(dzv.*v));
D_u=squeeze(sum(dzu));
D_v=squeeze(sum(dzv));
ubar=squeeze(hu./D_u);
vbar=squeeze(hv./D_v);
%
return
