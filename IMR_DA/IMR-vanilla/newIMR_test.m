% This script imports experimental data, for now just simuilated data with
% given parameters

addpath ./functions/
addpath ./src/

%% "True" parameters for simulation
timesteps = 90%10000; %90; % time steps 

%tspan = 2e-4;  % time to run simulation
tspan = 0.6e-4;
R0_true = 350.e-6; % mean soft %R0 = 100e-6; %R0 = 50.e-6; % Intital Radii
NT = 400; % Ammount of nodes inside the bubble (~100 is a good to start)
NTM = 20; % Ammount of nodes outside the bubble (~100 is a good to start)
%NT = 500;
%NTM = 500;
Pext_type = 'IC'; % Type of external forcing. Refer to RP_Cav
%Pext_Amp_Freq_true =[2000 0]; %Pext_Amp_Freq =[3550 0]; % [ Pressure ; Freq ]
disptime = 0; % 1 = display simulation time
Tgrad = 1; % Thermal effects inside bubble
Tmgrad = 1;% Thermal effects outside bubble
Cgrad = 1;  % Vapor diffusion effects
Dim = 0; % Output variables in dimensional form
comp = 1; % Activates the effect of compressibility
model = 'neoHook'; % 'neoHook', 'nhzen', 'sls' or 'linkv'

G1 = inf;
G_true = 7.69e3; % Estrada paper stiff
mu_true = 0.101; % Estrada paper stiff

% Added from old version:
alpha_true = 0;
lambda_nu_true = 0;
RelTol = 1e-7;
K=0; %for now

props.G = G_true;
props.mu = mu_true;
props.G1 = inf;
props.alpha = alpha_true;
props.lambda_nu = lambda_nu_true;

%G_true = 2.12e3; % Estrada paper soft
%mu_true = 0.118; % Estrada paper soft
Req_true = 0.15;

%G_true = 0; % water
%mu_true = 0.001; % water

Pmt_true = IMRcall_parameters(R0_true,G_true,G1,mu_true); % Calls parameters script
Ca_true = Pmt_true(5); Re_true = Pmt_true(6);
P_inf = Pmt_true(19); T_inf = Pmt_true(20);

ST = 0.056; % (N/m) Liquid Surface Tension
P_true = (P_inf + (2*ST)/(Req_true*R0_true) - Pvsat(T_inf))*(Req_true^3);

Pext_Amp_Freq_true = [P_true 0]; %as close to exp as possible

[t_true, R_true, Rdot_true, P_true, S_true, T_true, C_true, Tm_true, ...
    tdel_true, Tdel_true, Cdel_true,Tau_true] = funIMRsolver(model, props, ...
    tspan, R0_true, NT, NTM, Pext_type, Pext_Amp_Freq_true , disptime, Tgrad, ...
    Tmgrad, Cgrad, Dim, comp, Req_true,RelTol); %,timesteps)

xth = [R_true,Rdot_true,P_true,S_true,Tau_true,C_true,Tm_true]';

%% Add Ca and Re to truth vector

xth = [xth;log(Ca_true)*ones(size(R_true'));log(Re_true)*ones(size(R_true'))];

%% Adding noise to true value (simulated experimental measurements)
%std = 0.02;
yth = R_true' + std*randn(size(R_true'));

n = length(t_true)-1;

%% Getting needed values from data now

R0 = yth(1)*R0_true;