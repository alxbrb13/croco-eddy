% keys_sources_AVISO.m
%
%   keys_sources sets user defined paths and user keys
%
% Paths:
%   - path_in: directory containing the input files;
%              (default is '..\Data\')
%   - path_out: directory for the output files;
%              (default is '..\Tracks\')
%   - nc_u: full name of the netcdf file with the zonal component of
%           velocity (ssu) and the time index (step)
%   - nc_v: full name of the netcdf file with the meridional component of
%           velocity (ssv) and the time index (step)
%   - nc_dim: full name of the netcdf file with the domain coordinates 
%            (longitude and latitude) and the velocity mask (land-points=0;
%             ocean-points=1)
%
% User option keys:
%   - type_detection: flag to choose the field use as streamlines
%           1 : using velocity fields
%           2 : using ssh
%           3 : using both velocity fields and ssh for eddy contour, 
%               and keep the contour of the last one if it existS.
%   - extended_diags: flag to have extra diags concerning eddies directly
%       computed (like ellipse features or vorticity for each eddy)
%           0 : not computed
%           1 : computed as the same time as eddy detection
%           2 : computed afterward
%   - streamlines and daystreamfunction: save streamlines at steps 
%       'daystreamfunction' and profils of streamlines scanned as well
%       (1:stepF by default )
%   - periodic: flag to activate options for East-West periodic
%               (e.g. global fields or idealized simulations) domains.
%               IMPORTANT: first and last columns of the domain must be
%                          identical for this to work properly!!!!!
%
%-------------------------
% IMPORTANT - Input file requirements:
%
% All the variables are read from netcdf file.
% The package requires 3 different input files:
% 1) nc_dim with variables x(j,i),y(j,i) and mask(y,x)
% 2) nc_u with variable ssu(t,j,i) in m/s, step(t)
% 3) nc_v with variable ssv(t,j,i) in m/s, step(t)
% you can have also the sea level in an other file
% 4) nc_ssh with variable ssh(t,j,i) in m, step(t)
%
% t is the temporal dimension (number of time steps)
% j is the zonal dimension (number of grid points along y or latitude)
% i is the meridional dimension (number of grid points along x or longitude)
%
% The grid is assumed to be rectangular, with orientation N-S and E-W. 
% Grid indices correspond to geography, so that point (1,1) represents the
%      south-western corner of the domain.
% Latitudinal and longitudinal grid spacing can vary within the grid domain.
%
%-------------------------
%   Ver Sept 2021 Barboni Alexandre for CROCO output, integrate geostrophy
%   Ver Jun 2018 Briac and Romain Pennel
%   Ver Apr 2015 Briac Le Vu
%   Ver. 2.1 Oct.2012
%   Ver. 2.0 Jan.2012
%   Ver. 1.3 Apr.2011
%   Ver. 1.2 May.2010
%   Ver. 1.1 Dec.2009
%   Authors: Francesco Nencioli, francesco.nencioli@univ-amu.fr
%            Charles Dong, cdong@atmos.ucla.edu
%-------------------------
%
%=========================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User modification ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Experiment settings
%----------------------------------------------

% type of input sources (suffixe 1)
source = 'WindEddy';

% name for the experiment (suffixe 2)
config = '';

% use to diferenciate source field of surface height (adt, ssh, psi,...)
sshtype='croco'; % adt or sla

% use to differenciate job (prefixe 1)
runname ='QWA_rev1h';
expname='qwa_rev1h-2';
gridname='grd_lat385_1km.nc';

postname=runname;
% set the paths
path_in = ['/home6/datawork/abarboni/CROCO/RESULT/',runname,'/'];
path_out = ['/home6/datawork/abarboni/CROCO/RESULT/',runname,'/'];
%path_rossby = '/home6/datahome/abarboni/AMEDA/Rossby_radius/';

disp(['Compute from ',path_in])
disp([' to ',path_out])

if exist(path_out,'file')==0
    system(['mkdir ',path_out]);
    system(['mkdir ',path_out,'/']);
end

% add path_out path
addpath(path_out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input data file absolute name
nc_dim = [path_in,gridname];
nc_u = [path_in,expname,'_ssu.nc']; %Zonal speed file
nc_v = [path_in,expname,'_ssv.nc']; %Meridional speed file
nc_ssh = [path_in,expname,'_ssh.nc']; %ssh file

% Rossby deformation radius file
mat_Rd = [path_rossby,'']; % Not important, for regional simultion
%name_Rd = 'Rd_baroc1_extra';

% searched eddies typical radius
% eddies smaller than 1/4 this radius will be smoothed
Rd_typ = 11; % km

% minimal size for rmax to be reasonably detected (in number of Dx)
nRmin = 4; % 4 for 500m, 2 for 1km, A.Barboni changed : At 1/2 for AVISO

% In mod_eddy_centers :  length(radius)>=2 && radius(end)>=nRmin*Dx && radius(end)<=nR_lim*Rd

% variable names (could be automatised)
y_name = 'y_rho';
x_name = 'x_rho';
m_name = 'mask';
u_name = 'ugeo';%'u'
v_name = 'vgeo';%'v'
s_name = 'zeta';

% duration experiment (should read a 'time' variable)
if ~exist('stepF','var')
    u0 = squeeze(ncread(nc_ssh,s_name,[1 1 1],[1 1 Inf]));
    stepF = length(u0);
    clear u0
end

disp([' ',num2str(stepF),' time steps'])
disp(' ')

% rotation period (T) per day and time step in days (dps)
T = 3600*24; % day period in seconds
dps = 2; % number of day per time step   <===== CHANGE HERE IF RECORD TIMESTEP IS CHANGED

f_0=9.05e-5;  %  <===== CHANGE HERE IF LATITUDE IS CHANGED

level = 1; % vertical level (used only if any)

% degradation factor to test the algorithm
if ~exist('deg','var')
    deg = 2; % from 1 (default) to >10 in some experiment
end

%% Experiment option keys
%----------------------------------------------

% grid type
grid_ll = 0;
        % 0 : spatial grid in cartesian coordinates (x,y)
        % 1 : spatial grid in earth coordinates (lon,lat)

% grid regular or not (like arakawa in NEMO)
grid_reg = 1;
        % 0 irregular 
        % 1 regular 

% choose the field use as streamlines

type_detection = 1;
        % 1 : using velocity fields
        % 2 : using ssh
        % 3 : using both velocity fields and ssh, 
        %     and keep max velocity along the eddy contour

% if you want extended diags directly computed
extended_diags = 0;
        % 0 : not computed
        % 1 : computed as the same time as eddy detection
        % 2 : computed afterward  

% save streamlines at days daystreamfunction and profil as well
streamlines = 1;
daystreamfunction = 1:stepF;

% in case of periodic grid along x boundaries
periodic = 1;

%## New A. Barboni for periodic CROCO fields
borders=40 ; %% If periodic, number of PIXELS added at the boundary for enlarged fields
%Careful : in native grid pixel size, can be reduced by "deg" factor

% to keep firts and last detection after the tracking in NRT configuration
nrt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of user modification ---------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
