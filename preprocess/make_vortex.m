%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a CROCO VORTEX configuration.
%  Create a grid file and an intial file for the vortex experiment
%  (for the parent and the child grid).
%  The vortex is defined such as there is no motion below a definied
%  depth H0.
% 
%  Further Information:  
%  http://www.croco-ocean.org
%  
%  This file is part of CROCOTOOLS
%
%  CROCOTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  CROCOTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Ref: Penven, P., L. Debreu, P. Marchesiello et J.C. McWilliams,
%       Application of the ROMS embedding procedure for the Central 
%      California Upwelling System,  Ocean Modelling, 2006.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
addpath(genpath('/home6/datahome/abarboni/CROCO/croco-eddy/preprocess/')) %adding  with subfolders
cd('/home6/datahome/abarboni/CROCO/croco-eddy/preprocess/')

%  Title
%
title='VORTEX';
%
%  Names
%
parent_grd='vortex_grd.nc';
parent_ini ='vortex_ini.nc';
%
% Parameters for the Parent Grid
%
dx=2e3;                % Horizontal resolution  %%% CHANGED from initial resol 
xmax=100e3;             % Domain length
H0=3000;                % Depth
H=50;                 % Level of no-motion
theta=38.5; %38.5            % Latitude (beta-plane)
R=6367442.76;           % Earth radius
Pa=1013e2;              % Atmospheric pressure
rho0=1029;            % Mean ocean density
umax=0.4;                 % Max velocity (<0 cyclonic in the northern hemisphere)
%lambda=sqrt(2)*25e3;    % Vortex radius (exp(r2/lambda2))
Rmax=25e3;
alpha=1.6;
g=9.81;                 % Gravity acceleration
%N2=(0.003)^2;           % Brunt-Vaissala frequency

%%% Strati params for rho
R1=29.03;Rho_0=25.3;Sl=1.8e-4;ZT=55;
%%%%%% Strati params for temperature
%T1=16.8; Tsurf=30.3;ZT=45; St=3.5e-3

% Perturbation params
ModeMax=10;              % Max radial wavenumber introduced 
Perturb=0;               % Perturb factor (fraction of initial tangential speed)
% keep Perturb=0 unless you want to test stability analysis
%CycleVert=0;       % number of vertical nodes  introduced (Not active)

% Vertical grid parameters
%
N=80;  
theta_s=8;
theta_b=0;
hc=400;
vtransform =  2.; % s-coordinate type (1: old- ; 2: new- coordinates):wq:

%
% Nesting parameters
%
refinecoeff=3;
jmin=21;
jmax=40;
imin=21;
imax=40; 
%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
%
% Horzontal Grid
%
x=[-xmax-dx/2:dx:xmax+dx/2];
y=x;
dy=dx;
[X,Y]=meshgrid(x,y);
%
% Topo
%
h0=H0+0*X;
%
% Coriolis term (beta plane)
%
deg2rad=pi/180;
omega=2*pi/(24*3600);
f0=2*omega*sin(deg2rad*theta);
beta=2*omega/R*cos(deg2rad*theta);
f=f0; %+beta*Y;
%
% Compute zeta,ubar,vbar,u,v,t for the vortex
%
barocvortex_modif
%
%  Create the grid file
%
[Mp,Lp]=size(zeta);
M=Mp-1;
L=Lp-1;
nc=netcdf(parent_grd, 'clobber');
%%redef(nc);            % for Octave compatiblity
nc('xi_rho') = Lp;
nc('eta_rho') = Mp;
nc('xi_psi') = L;
nc('eta_psi') = M;
nc('one') = 1;
nc{'el'} = ncdouble('one');
nc{'xl'} = ncdouble('one');
nc{'spherical'} = ncchar('one');
nc{'h'} = ncdouble('eta_rho', 'xi_rho');
nc{'f'} = ncdouble('eta_rho', 'xi_rho');
nc{'pm'} = ncdouble('eta_rho', 'xi_rho');
nc{'pn'} = ncdouble('eta_rho', 'xi_rho');
nc{'x_rho'} = ncdouble('eta_rho', 'xi_rho');
nc{'y_rho'} = ncdouble('eta_rho', 'xi_rho');
nc{'mask_rho'} = ncdouble('eta_rho', 'xi_rho');
%%endef(nc);		% for Octave compatiblity
nc.title = ncchar(title);
nc.title = title;
nc.date = ncchar(date);
nc.date = date;
nc.type = ncchar('CROCO grid file');
nc.type = 'CROCO grid file';
%
%  fill the grid file
%
nc{'xl'}(:)=dx*(L-1);
nc{'el'}(:)=dy*(M-1);
nc{'spherical'}(:)='F';
nc{'h'}(:)=h0;
nc{'f'}(:)=f;
nc{'pm'}(:)=1/dx;
nc{'pn'}(:)=1/dx;
nc{'x_rho'}(:)=X;
nc{'y_rho'}(:)=Y;
nc{'mask_rho'}(:)=1+0.*Y;
close(nc);
%
%  Create and fill the initial file
%
create_inifile(parent_ini,parent_grd,title,...
               theta_s,theta_b,hc,N,0,'clobber',vtransform)
nc=netcdf(parent_ini,'write');
nc{'u'}(:) =  u; 
nc{'v'}(:) =  v; 
nc{'zeta'}(:) =  zeta; 
nc{'ubar'}(:) =  ubar; 
nc{'vbar'}(:) =  vbar; 
nc{'temp'}(:) =  t; 
close(nc)


