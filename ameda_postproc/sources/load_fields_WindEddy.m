function [x,y,mask,u,v,ssh] = load_fields_CROCO(stp, resolution, degra)
%[x,y,mask,u,v,ssh]=load_field(stp {,resolution,degra});
%
%  Get netcdf and load the grid and the velocities field (also ssh if any)
%  - degrad ('deg') is a sample factor (use in some experiment).
%  - resolution ('resol') is the factor to interpolate the grid by a resolution factor res.
%   res=0: interpolation on regular grid if necesseray
%   res>1: interpolate the grid by a resolution factor res
%
%  Enlarge the mask into land by adding b pixels of nil velocities into the
%  land and 1 pixel in land from an averaged of 9 neighbours for ssh
%	u(b)=0,u(~mask(b))=nan;
%	v(b)=0,v(~mask(b))=nan;
%	ssh(b)=mean(ssh)(b-1) or mean(ssh),ssh(~mask(b-1))=nan;
%
%  Output are matlab matrice used with tracking_plot routines and should be
%  saved in fields.mat file and in fields_inter.mat.
%  fields size must be [lat,lon,time]
%  output velocities must be in m/s and ssh in m
%
%  For a description of the input parameters see mod_eddy_param.m.
%
%-------------------------
%  June 2018 Briac Le Vu and Romain Pennel
%  Sept 2021 Alexandre Barboni => Periodic N-S + W-E boundaries, basically
%       pass through initial field u0 => Periodic up => degraded (or not) u
%       => interpolated ui. Assumes u/v grid is identical (no rho/psi diff grid)
%-------------------------
%
%=========================

% load keys_sources and parameters (use mod_eddy_params.m first)
%%%%%%%%%%%%%%%%%%%%%%%% Begin Modif ABarboni 
%borders=5%
run(['keys_sources_WindEddy_'])  % Not necessary, just to load param borders
%, because not in param_eddy_tracking
%%%%%%%%%%%%%%%%%%%%%%%% End Modif ABarboni 
load('param_eddy_tracking')


%----------------------------------------
% replace parameters by arguments
if nargin==3
% No interpolation by default
    deg = degra;
    resol = resolution;
elseif nargin==2
% No degradation and No interpolation by default
    resol = resolution;
end

%% load fields
disp(['Get grid and velocities field at step ',num2str(stp),' ...'])

lon0 = double(ncread(nc_dim,x_name))';
lat0 = double(ncread(nc_dim,y_name))';
%elseif strcmp(source,'CROCO')   %%% MODIF ABarboni
mask0 = ones(size(lat0))';  
%mask0 = double(ncread(nc_dim,m_name))'; %%% End Modif ABarboni

%----------------------------------------
% interpolate u and v
u0 = double(squeeze(permute(ncread(nc_u,u_name,[1 1 stp],[Inf Inf 1]),[2,1,3])));
v0 = double(squeeze(permute(ncread(nc_v,v_name,[1 1 stp],[Inf Inf 1]),[2,1,3])));

if type_detection>=2
    ssh0 = squeeze(permute(ncread(nc_ssh,s_name,[1,1,stp],[Inf,Inf,1]),[2,1,3]));
else
    ssh = [];
end
%----------------------------------------
lon0(isnan(lon0))=0;
lat0(isnan(lat0))=0;
u0(isnan(u0))=0;
v0(isnan(v0))=0;

%%%%%%%%%%%%%%%%%%%%%%%% Begin Modif ABarboni 
[N0,M0] = size(u0);

%% Enlarge degraded field % % New ABarboni
if periodic
    % add box pixels to the boundaries (East-West periodicity) AND
    % North-South periodicity
    % u and v matrices are expanded by adding the ending columns
    % to the beginning of the domain, and the beginning columns at its end.
    up=zeros(N0+2*borders, M0+2*borders);
    vp=zeros(N0+2*borders, M0+2*borders);
    latp=zeros(N0+2*borders, M0+2*borders);
    lonp=zeros(N0+2*borders, M0+2*borders);
    maskp=zeros(N0+2*borders, M0+2*borders);
    
    % enlarged u,v with the same values at WEST and EAST boundary
    up(borders+1:N0+borders,:)=[u0(:,end-borders:end-1),u0,u0(:,2:borders+1)];
    vp(borders+1:N0+borders,:)=[v0(:,end-borders:end-1),v0,v0(:,2:borders+1)];
    maskp(borders+1:N0+borders,:)=[mask0(:,end-borders:end-1),mask0,mask0(:,2:borders+1)];
    
    Dlon=lon0(1,1)-lon0(1,end);
    Dlat=lat0(1,1)-lat0(end,1);
    lonp(borders+1:N0+borders,:)=[lon0(:,end-borders:end-1)+Dlon,lon0,lon0(:,2:borders+1)-Dlon];
    latp(borders+1:N0+borders,:)=[lat0(:,end-borders:end-1),lat0,lat0(:,2:borders+1)];
    
    % enlaregd u,v with the same values at NORTH and SOUTH boundary
    up(1:borders,:)=up(N0:end-borders-1,:);%up(N0+1:end-borders,:);
    up(N0+borders+1:end,:)=up(borders+2:2*borders+1,:);%up(borders+1:2*borders,:);
    vp(1:borders,:)=vp(N0+1:end-borders,:);
    vp(N0+borders+1:end,:)=vp(borders+1:2*borders,:);
    
    latp(1:borders,:)=latp(N0:end-borders-1,:)+Dlat;
    latp(N0+borders+1:end,:)=latp(borders+2:2*borders+1,:)-Dlat;
    
    lonp(1:borders,:)=lonp(N0:end-borders-1,:);
    lonp(N0+borders+1:end,:)=lonp(borders+1:2*borders,:);
    
    maskp(1:borders,:)=maskp(N0:end-borders-1,:);
    maskp(N0+borders+1:end,:)=maskp(borders+1:2*borders,:);
else
    % add no pixels otherwise
    borders = 0;
    up=u0;
    vp=v0;
    latp=lat0;
    lonp=lon0;
    maskp=mask0;
end
%%%%%%%%%%%%%%%%%%%%%%%% End Modif ABarboni
%% Work on degraded field 
if deg~=1
    disp(['  Fields are degraded by a factor ',num2str(deg)])
    disp('  (degraded grid becomes native grid)')
end

% produce degraded field 
%%%%%%%%%%%%%%%%%%%%% Begin Modif A. Barboni
if strcmp(source,'WindEddy')   
    x = lonp(1:deg:end,1:deg:end)/1000;
    y = latp(1:deg:end,1:deg:end)/1000;
else
    x = lonp(1:deg:end,1:deg:end);
    y = latp(1:deg:end,1:deg:end);
end
mask = maskp(1:deg:end,1:deg:end);
u = up(1:deg:end,1:deg:end);
v = vp(1:deg:end,1:deg:end);
if type_detection>=2
    ssh = sshp(1:deg:end,1:deg:end);
end
%%%%%%%%%%%%%%%%%%%%%% END Modif A.Barboni

% get the grid size
[N,M] = size(x);


%% Increase resolution r factor by linear interpolation
if resol==1 && grid_reg

    disp('NO INTERPOLATION')
    
    % fix fields to NaN in land 
    u(mask==0) = NaN;
    v(mask==0) = NaN;

    % Enlarge mask into land by 1 pixel and compute ssh in the first land pixel if needed
    disp('Enlarge coastal mask by adding 1 pixel of ocean to the coast ...')
    for i=1:N
        for j=1:M
            if mask(i,j)==0 &&...
                    sum(sum(mask(max(i-1,1):min(i+1,N),max(j-1,1):min(j+1,M))))~=0
                u(i,j) = 0;
                v(i,j) = 0;
                if type_detection>=2 && isnan(ssh(i,j))
                    ssh1 = ssh(max(i-1,1):min(i+1,N),max(j-1,1):min(j+1,M));
                    ssh(i,j) = nanmean(ssh1(:));
                end
            end
        end
    end
    
elseif ( resol==1 && ~grid_reg ) || resol ~= 1

    if grid_reg
        disp(['Change resolution by computing 2D SPLINE INTERPOLATION ',...
        'by a factor ',num2str(resol)])
    else
        disp('No change in resolution, REGRIDDING from Arakawa to regular grid')
    end
    
    %----------------------------------------
    % Increase resolution of the mask
    [Ni,Mi] = size(xi);
    if size(x,1) ~= size(u,1) && Ni ~= resol*(N-1)+1
        disp(['Need to adapt gridvel.mat to deg=',num2str(deg),' and res=',num2str(resol)])
        return
    end
    
    %----------------------------------------
    % Enlarge mask into land by 1 pixel and compute ssh in the first land pixel if needed
    maski1 = maski;% enlarged mask
    disp('Enlarge coastal mask by adding 1 pixel of ocean to the coast ...')
    for i=1:Ni
        for j=1:Mi
            if maski(i,j)==0 &&...
                    sum(sum(maski(max(i-1,1):min(i+1,Ni),max(j-1,1):min(j+1,Mi))))~=0
                maski1(i,j)=1;
            end
        end
    end

    %----------------------------------------
    % fix fields to 0 in land 
    u(mask==0 | isnan(u)) = 0;
    v(mask==0 | isnan(v)) = 0;
    % Increase resolution of fields (griddata with irregular grid)
    ui = griddata(x,y,u,xi,yi,'cubic');
    vi = griddata(x,y,v,xi,yi,'cubic');
    if type_detection>=2
        sshi = griddata(x,y,ssh,xi,yi,'cubic');
    else
        sshi = [];
    end
    % Mask velocities and ssh with their enlarged mask
    ui(maski1==0) = NaN;
    vi(maski1==0) = NaN;
    if type_detection>=2
        sshi(maski1==0) = NaN;
    end

    %----------------------------------------
    % Export interpolated fields
    x = xi;
    y = yi;
    u = ui;
    v = vi;
    mask = maski;
    ssh = sshi;
       
end
