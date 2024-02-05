function [ugeo,vgeo] = compute_speed_geo_from_grd(n,f,x,y,ssh)
% function [ugeo,vgeo] = compute_speed_geo(n,x,y,ssh)
%
%  compute geostrophique field by a n-stencil centered finite difference 
%  with non-uniform grid spacing in x and y from a ssh field
%
% September 2021 Alexandre Barboni to compute from cartesian grid
%-------------------------
%  April 2020 Briac Le Vu nased on https://github.com/mluhar/dynamicblade 
%-------------------------
%
%=========================

%----------------------------------------

% gravitation constant
g = 9.8; % m.s-2

%----------------------------------------
% Calculation of finite spatial element
[~,dx] = gradient(x);
[dy,~] = gradient(y);

% Calcul finite spatial element in km
%dcy=cumsum(dy,2)*R;
%dcx=cumsum(dx*R.*cosd(y),1);

% in meters
dcy = cumsum(dy,2);%dcx*1000; % m
dcx = cumsum(dx,1);%dcy*1000; % m

%----------------------------------------
% initialisation
ugeo = nan(size(ssh));
vgeo = nan(size(ssh));

%----------------------------------------
disp('Compute geostrophic speed')
for t=1:size(ssh,3)

  disp(['  step ',num2str(t),'/',num2str(size(ssh,3))])

  % u from n-stencil centered difference of dy
  for i=1:size(dcy,1)
  
    % mask land and island along y(i,:)
    idnan = isnan(ssh(i,:,t));

    if length(dcy(i,~idnan))>n

      % the n dy Weight at every y(i,j)
      Dynan = fdmatrix(dcy(i,~idnan),1,n);
  
      % number of remaining ssh(i,:)
      Nnan = length(ssh(i,~idnan,t));

      % prepare the ssh matrix of ssh along y(i,:)
      sshi = repmat(ssh(i,~idnan,t),[Nnan 1]);

      % u(i,j) is the sum of Dy weighted ssh along y(i,:)
      %ugeo(i,~idnan,t) = -g./f(i,~idnan).*nansum(Dynan.*sshi,2)';
      ugeo(i,~idnan,t) = -g./f.*nansum(Dynan.*sshi,2)';
    end

  end

  % v from n-stencil centered difference of dx latitude by latitude j
  for j=1:size(dcx,2)
  
    % mask land and island along x(:,j)
    idnan = isnan(ssh(:,j,t));

    if length(dcx(~idnan,j))>n

      % the n dx Weight at every x(i,j)
      Dxnan = fdmatrix(dcx(~idnan,j),1,n);
  
      % number of remaining ssh(:,j)
      Nnan = length(ssh(~idnan,j,t));
  
      % prepare the ssh matrix of ssh along x(:,j)
      sshj = repmat(ssh(~idnan,j,t)',[Nnan 1]);
  
      % v(i,j) is the sum of Dx weighted ssh along x(:,j)
      %vgeo(~idnan,j,t) = g./f(~idnan,j).*nansum(Dxnan.*sshj,2);
      vgeo(~idnan,j,t) = g./f.*nansum(Dxnan.*sshj,2);

    end

  end

end

end % end main function




