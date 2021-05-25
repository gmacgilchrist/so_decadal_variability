clear; close all;
p=genpath('~/packages/gcmfaces'); addpath(p);
eccodir = '/data2/project/ECCO4v4/'
if ~exist('nctiles_grid','dir')
	command = ['ln -s ' eccodir 'nctiles_grid nctiles_grid'];
	status = system(command);
end
global mygrid
grid_load
gcmfaces_global;

yr = 1992; mt = 01;
suffix = [num2str(yr) '_' num2str(mt,'%02i') '.nc'];
myenv.nctilesdir = fullfile(eccodir,['GAMMAN/' num2str(yr) '/'])
fld=read_nctiles([myenv.nctilesdir 'GAMMAN']);

%myenv.nctilesdir = fullfile(eccodir,['nctiles_monthly/GAMMAN/' num2str(yr) '/'])
%fld=read_nctiles([myenv.nctilesdir 'GAMMAN']);
%GAMMAN = ncread([eccodir 'GAMMAN/' num2str(yr) '/GAMMAN_' suffix],'GAMMAN');
%THETA = ncread([eccodir 'nctiles_monthly/THETA/' num2str(yr) '/THETA_' suffix],'THETA');
%SALT = ncread([eccodir 'nctiles_monthly/SALT/' num2str(yr) '/SALT_' suffix],'SALT');

%dTdx = calc_T_grad(THETA)
 
