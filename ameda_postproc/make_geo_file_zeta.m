%%%%%
%
%  Script to compute SSH from CROCO output, modified from Briac Le Vu
%  script, september 2021
%
%%%%%
start
clear all
close all

addpath('/home6/datahome/abarboni/AMEDA/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User modifications
expname='qwa_rev1h-2'
filename=[expname,'_ssh.nc'];
gridname='grd_lat385_1km.nc'

% directory
dirin=['/home6/datawork/abarboni/CROCO/RESULT/QWA_rev1h/'];

n=9; % n-stencil finite difference centred method
% Coriolis parameter
f=9.05e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
% read fields
ssh=ncread([dirin,filename],'zeta');
[N,M,L]=size(ssh);

%lon=ncread([dirin,'lon_lat_',sshname,'_',domname,'.nc'],'lon');
%lat=ncread([dirin,'lon_lat_',sshname,'_',domname,'.nc'],'lat');

lon=ncread([dirin,gridname],'x_rho');
lat=ncread([dirin,gridname],'y_rho');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute geostrophic velocities stencil-9, script with cartesian grid

[ugeo,vgeo] = compute_speed_geo_from_grd(n,f,lon,lat,ssh);

nccreate([dirin,expname,'_ssu.nc'],'day','Dimension',{'day' Inf});
nccreate([dirin,expname,'_ssu.nc'],'ugeo','Dimension',{'x' N 'y' M 'day' L});
ncwrite([dirin,expname,'_ssu.nc'],'ugeo',ugeo);

nccreate([dirin,expname,'_ssv.nc'],'day','Dimension',{'day' Inf});
nccreate([dirin,expname,'_ssv.nc'],'vgeo','Dimension',{'x' N 'y' M 'day' L});
ncwrite([dirin,expname,'_ssv.nc'],'vgeo',vgeo);





