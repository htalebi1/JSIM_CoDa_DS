%----------------------------------------------------------------------------------------------------------------
% Authors: Hassan Talebi, Ute Mueller, Raimon Tolosana-Delgado
% Paper: "Joint simulation of compositional and categorical data via direct sampling technique – Application to
%         improve mineral resource confidence"
% Journal: Computer & Geosciences
%----------------------------------------------------------------------------------------------------------------

A MATLAB code (DS.m) has been developed for the joint simulation of compositional and categorical variables. The following 
parameters together with a partially informed grid should be used as inputs to this code. 

%%%--------------------------
% Simulation parameters
%%%--------------------------
n_cat = 1;                                         % number of categorical variables
n_cont = 4;                                        % number of continuous variables
n = 12;                                            % maximum number of points in the data event
f = 0.8;                                           % maximum fraction of the training image to be scanned
t = 0.025;                                         % overall distance threshold (between 0 and 1)
n_realiz = 1;                                      % number of realizations
postpone=1;                                        % times to postpone the simulation of a node when a close pattern not found
nx=250;                                            % number of nodes in x direction
ny=250;                                            % number of nodes in y direction
nz=1;                                              % number of nodes in z direction
dimx=1;                                            % cell dimension in x direction
dimy=1;                                            % cell dimension in y direction
dimz=1;                                            % cell dimension in z direction
expo_w=2;                                          % exponent of the inverse distance weights
cond_w=5;                                          % weight of hard data vs simulated data
alpha_w=0.5;                                       % mixing coefficient 
var_w=[0.25,0.25,0.25,0.25];                       % weight of each variable (first categorical variables 
rx=50;                                             % maximum search distances along x direction
ry=50;                                             % maximum search distances along y direction
rz=5;                                              % maximum search distances along z direction
denoise=1;                                         % noise removal (resimulating via fully informed nodes)
ilr_trans = true;                                  % transform compositional data to ilr space
plt=true;                                          % plot