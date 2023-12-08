% This file sets up an initial ensemble for experimental data setups
% State vector used: [R,U,P,S,Tau,C,Tm,log(Ca),log(Re)]
addpath ./IMR-vanilla/functions
addpath ./IMR-vanilla/src
% Shuffle rng to ensure randomness of results
rng('shuffle')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Set EnKF parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Guess for parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Req and calculate initial partial pressure
%ST = 0.056; % (N/m) Liquid Surface Tension
Pmt_temp = IMRcall_parameters(R0,G_guess,G1_guess,mu_guess); % Calls parameters script
Ca = Pmt_temp(5); Re = Pmt_temp(6);
P_inf = Pmt_temp(19); T_inf = Pmt_temp(20);
Req = mean(yth(end-5:end)); % take mean of end of sim to be Req
P_guess = (P_inf + (2*ST)/(Req*R0) - Pvsat(T_inf))*(Req^3);
Pext_Amp_Freq =[P_guess 0]; % [ Pressure ; Freq ], Freq = 0 for IC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Determine initial state vector based on parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialize
x_init = x0_true';
spread = [Rspread; Uspread; Pspread; Sspread; ones(NT,1)*tauspread; ...
    ones(NT,1)*Cspread; ones(NTM,1)*Tmspread; Brspread; fohspread; ...
    Caspread; Respread; Despread; alphaspread; lambda_nuspread];
xi = (1 + spread .* randn(N,q)) .* repmat(x_init,1,q) + ...
    repmat([0;0;0;0;zeros(2*NT+NTM,1);0;0;0;0;0;0;0],1,q) .* randn(N,q);
xi(3,:) = log(xi(3,:));
xi(3,:) = log(x0_true(3));
%xi = [xi(1:2*NT+NTM+4,:);log(xi(2*NT+NTM+5:end,:))];
xi = [xi(1:2*NT+NTM+6,:);log(xi(2*NT+NTM+7:2*NT+NTM+8,:));xi(end-2,:);log(xi(end-1,:));xi(end,:)]; % for now
%%% Other parameters to initialize
x_est = zeros(N,n+1);