% This file sets up an initial ensemble for experimental data setups

% State vector used: [R,U,P,S,Tau,C,Tm,log(Ca),log(Re)]

addpath ./IMR-vanilla/functions
addpath ./IMR-vanilla/src

% Shuffle rng to ensure randomness of results
rng('shuffle')

%% Set EnKF parameters
% Ensemble size
%q = 10;

% Covariance Inflation parameters
%beta = 1.02; % 1.005 <= beta <= 1.05
%alpha = 0.005; % *random noise in forecast step
%CI_theta = 0.6; % 0.5 <= theta <= 0.95
%CI_scheme = 2; % 1 is scalar alpha, 2 is relaxation prior to mean
%CI_add = 1; % set to 1 for additive covariance

%init_scheme = 2; % 1 is POD method, 2 is truth+noise IN OLD VERSION

%% Guess for parameters

%NT = 500; %400 % Ammount of nodes inside the bubble (>=500 is a good to start)
%NTM = 500; %20 % Ammount of nodes outside the bubble (>=10 is a good to start)
%Pext_type = 'IC'; % Type of external forcing. Refer to RP_Cav

% Find Req and calculate initial partial pressure
% Needed parameters
%ST = 0.056; % (N/m) Liquid Surface Tension
Pmt_temp = IMRcall_parameters(R0,G_guess,G1_guess,mu_guess); % Calls parameters script
Ca = Pmt_temp(5); Re = Pmt_temp(6);
P_inf = Pmt_temp(19); T_inf = Pmt_temp(20);

Req = mean(yth(end-5:end)); % take mean of end of sim to be Req
P_guess = (P_inf + (2*ST)/(Req*R0) - Pvsat(T_inf))*(Req^3);

Pext_Amp_Freq =[P_guess 0]; % [ Pressure ; Freq ], Freq = 0 for IC
%Pext_Amp_Freq = [100 0];

% Simulation parameters
%disptime = 0; % 1 = display simulation time
%Tgrad = 1; % Thermal effects inside bubble

%if exist('Tmgrad')
%    % don't modify if already set
%else
%    Tmgrad = 1;% Thermal effects outside bubble
%end

%Cgrad = 1;  % Vapor diffusion effects
%Dim = 0; % Output variables in dimensional form
%comp = 1; % Activates the effect of compressibility

%
%G1 = inf;
%G = G_guess;
%mu = mu_guess;

%% Determine initial state vector based on parameters
initialize

%% Create initial ensemble
%{
% Use a certain spread around initial state vector given by initialize
x_init = (1 + 0.1*randn(N,q)) .* repmat(x0_true',1,q) + ...
    repmat([0;0;0;0;zeros(2*NT+NTM,1);50;70],1,q) .* randn(N,q);

xi = [x_init(1:end-2,:);log(x_init(end-1,:));log(x_init(end,:))];
xi(3,:) = log(xi(3,:));

% Maybe:
%xi(1,:) = 1;
xi(3,:) = log(x0_true(3));
%xi = [x_init(1:end-2,:);-Inf(1,q);log(x_init(end,:))];
%}

% Spread after transform:
%{
x_init = x0_true';
x_init(3) = log(x_init(3));
x_init = [x_init(1:end-2);log(x_init(end-1));log(x_init(end))];

xi = (1 + 0.001*randn(N,q)) .* repmat(x_init,1,q) + ...
    repmat([0;0;0;0;zeros(2*NT+NTM,1);0;2],1,q) .* randn(N,q);
%}

% Custom spread:

x_init = x0_true';

% New version:
%{
Rspread = 0.1;
Uspread = 0.1;
Pspread = 0.1;
Sspread = 0.1;
tauspread = 0.1;
Cspread = 0.001;
Tmspread = 0.0005;
Brspread = 0.01;
fohspread = 0.01;
Caspread = 0.25;
Respread = 0.25;
Despread = 0;
alphaspread = 0;
lambda_nuspread = 0;
%}
% Testing with (almost) no spread:
%{
Rspread = 0;
Uspread = 0;
Pspread = 0;
Sspread = 0;
tauspread = 0;
Cspread = 0.;
Tmspread = 0;
Brspread = 0;
fohspread = 0;
Caspread = 0;
Respread = 0;
Despread = 0;
alphaspread = 0;
lambda_nuspread = 0;
%}
%{
Rspread = 0.1;
Uspread = 0.1;
Pspread = 0.1;
Sspread = 0.1;
tauspread = 0.1;
Cspread = 0.01;
Tmspread = 0.1;
Caspread = 0.25;
Respread = 0.25;

% no spread test
%
Rspread = 0.1;
Uspread = 0.1;
Pspread = 0.1;
Sspread = 0.1;
tauspread = 0.1;
Cspread = 0.01;
Tmspread = 0.0001;
Tmspread = 0;
Caspread = 0.2;
Respread = 0.2;
%}

%Caspread = 0;
%Respread = 0;

%{
% other try for EnKF:
Rspread = 0.1;
Uspread = 0.1;
Pspread = 0.1;
Sspreat = 0.1;
tauspread = 0.01;
Cspread = 0.01;
Tmspread = 0.01;
Caspread = 0.1;
Respread = 0.1;
%}

spread = [Rspread; Uspread; Pspread; Sspread; ones(NT,1)*tauspread; ...
    ones(NT,1)*Cspread; ones(NTM,1)*Tmspread; Brspread; fohspread; ...
    Caspread; Respread; Despread; alphaspread; lambda_nuspread];

%xi = (1 + spread .* randn(N,q)) .* repmat(x_init,1,q) + ...
%    repmat([0;0.1;0.2;0.5;zeros(2*NT+NTM,1);0.01;0.01;0;0;0;0.01;0.01],1,q) .* randn(N,q);

xi = (1 + spread .* randn(N,q)) .* repmat(x_init,1,q) + ...
    repmat([0;0;0;0;zeros(2*NT+NTM,1);0;0;0;0;0;0;0],1,q) .* randn(N,q);

%xi = (1 + spread .* randn(N,q)) .* repmat(x_init,1,q) + ...
%    repmat([0;0.01;0.01;0.01;zeros(2*NT+NTM,1);0;0],1,q) .* randn(N,q);

%xi = (1 + spread .* randn(N,q)) .* repmat(x_init,1,q) + ...
%    repmat([0;0;0;0;zeros(2*NT+NTM,1);0;0],1,q) .* randn(N,q);

% Trying a temp modification to increase variance in ensemble members
% without crashing

%xi(2*NT+5:2*NT+NTM+4,:) = xi(2*NT+5:2*NT+NTM+4,:) + repmat(randn(1,q),NTM,1).*0.001;

xi(3,:) = log(xi(3,:));
xi(3,:) = log(x0_true(3));
%xi = [xi(1:2*NT+NTM+4,:);log(xi(2*NT+NTM+5:end,:))];
xi = [xi(1:2*NT+NTM+6,:);log(xi(2*NT+NTM+7:2*NT+NTM+8,:));xi(end-2,:);log(xi(end-1,:));xi(end,:)]; % for now
%

%% Other parameters to initialize

x_est = zeros(N,n+1);

