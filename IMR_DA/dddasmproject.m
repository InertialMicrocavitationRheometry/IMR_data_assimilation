% Data assimilation master file
% Author: Jean-Sebastien Spratt -- jspratt@caltech.edu

% This file is a wrapper which runs a variety of data assimilation methods
% on laser-induced cavitation radius vs time data. It uses the code from
% Spratt et al. (2020). Depending on which data assimilation method and
% physical model is used, parameters from different sections below must be
% specified
close all
clear all
clc

%% Data import
% see import_data_exp.m file for details on data loading. The following is
% meant for data in the same format as RofT data from Selda's experiments.
% If the format is different, the import_data_exp.m file will need to be
% modified accordingly

data_type = 'exp'; % 'sim' or 'exp'
data_set = 'SoftPA_nobeads';
data_filepath = (['example_data/']);
data_filename = 'normalized_unscaled_Experimental_Data_10PA_.06BIS.mat'; %'PolyAcry_12_48_08_updated.mat';% name of file containing R vs T data
num_peaks = 2; % number of 'peaks' to assimilate in radius data
               % (peaks = collapse points as in Estrada paper)

               
               
%% Data assimilation parameters

method = 'EnKS'; % data assimilation method ('En4D','EnKS',('EnKF'))

G_guess = 500;
G1_guess = 1e9;
mu_guess = 0.05;
alpha_guess = 0.5;
lambda_nu_guess = 0.1;
%}
q = 48; % Number of ensemble members
std = 0.01; % expected standard deviation of measurements;
init_scheme = 2; % leave as 2, initializes ensemble with truth + noise

% The following are ending criteria for iterative optimization:
epsilon = 1e-5; % threshold difference in norm of state vector between steps
max_iter = 5; % max # of iterations until end optimization (5 to 10 is usually good)

% IEnKS only:
if method == 'EnKS'
    l = 1   ; %lag of smoother
end

% Note: Beta coefficients are fixed as equal here. They can be modified in
% the main_En4D_peaks and main_mda_peaks files if needed, if more weight
% needs to be attributed to earlier or later points in assimilation

%% Modeling parameters
model = 'neoHook'; % 'neoHook','nhzen','sls','linkv','fung','fung2','fungexp','fungnlvis'
NT = 30; % Amount of nodes inside the bubble
NTM = 30; % Amount of nodes outside the bubble
Pext_type = 'IC'; % Type of external forcing
ST = 0.056; % (N/m) Liquid Surface Tension

Tgrad = 1; % Thermal effects inside bubble
Tmgrad = 0; % Thermal effects outside bubble
Cgrad = 1; % Vapor diffusion effects
comp = 1; % Activates the effect of compressibility (0=Rayleigh-Plesset, 1=Keller-Miksis)
disp_timesteps = 1; % displays timesteps in En4D run (for debugging)

% following should not be changed (untested):
disptime = 0; % 1 = display simulation time
Dim = 0; % 1 = output variables in dimensional form

%% Covariance Inflation parameters
% The following is only used for the IEnKS. Default values provided below
% See Spratt et al. (2020) section 3.3 for details on theta, lambda

CI_scheme = 2; % 1 is scalar alpha CI, 2 is RTPS (BEST)
CI_theta = 0.7; % multiplicative covariance parameter (0.5 < theta < 0.95)
CI_add = 0; % Set to 1 for additive covariance (else 0)
beta = 1.02; % additive covariance parameter (lambda in paper) (1.005 < beta < 1.05)
alpha = 0.005; % random noise in forecast step (only for EnKF)

%% Spread of parameters in the ensemble
%
%Rspread = 0.02;
Rspread = 0.01;
Uspread = 0.1;
Pspread = 0.1;
Sspread = 0.1;
tauspread = 0.1;
Cspread = 0.001;
Tmspread = 0.0005;
%Tmspread = 0.001;
Brspread = 0.01;
fohspread = 0.01;
%Caspread = 0.3;
Caspread = 0;
Respread = 0.3;
Despread = 0; % set to 0 if not used in model
alphaspread = 0.3; % set to 0 if not used in model
lambda_nuspread = 0; % set to 0 if not used in model
%}


%% Do not modify
visco_params = struct('G',G_guess,'G1',G1_guess,'mu',mu_guess, ...
    'alpha',alpha_guess,'lambda_nu',lambda_nu_guess);
est_params = [];
R = std^2; % Measurement error covariance

%% Run main for corresponding method:

if method == 'En4D'
    main_En4D_peaks
elseif method == 'EnKS'
    main_mda_peaks
end
%
% Save results
save(['./results/',data_filename(1:8),'_',method, ...
    '_',model,'_dataset',num2str(dataset),'_q',num2str(q),'_G', ...
    num2str(G_guess),'_mu',num2str(mu_guess),'_',num2str(num_peaks), ...
    'peaks',datestr(now),'.mat'],'-v7.3')
%}
%% Plotting

if method == 'En4D'
    plot_En4D_exp
elseif method == 'EnKS'
    plot_Kalman_exp
end

