% based on Sakov 2011 paper
% Author: Jean-Sebastien Spratt (jspratt@caltech.edu)

% IEnKS-MDA script

%rng(99)
tic
%% Clean up
%close all
%clear all
clc

%% Initialize and import data

% Some variables to set
%q = 10;
%timesteps = 301;

%import_data
if strcmp(data_type,'sim') == 1
    import_data_sim
elseif strcmp(data_type,'sim_file') == 1
    load(data_filename)
else
    import_data_exp
end

% Guess for G and mu: 
%G_guess = 6000;
%mu_guess = 0.08;

% Create Initial Ensemble
%create_ensemble;
create_ensemble_exp

tau_del = cell(q,1);
%% Iterate

vars = {NT Pext_type Pext_Amp_Freq disptime Tgrad Tmgrad ...
    Cgrad comp t0 neoHook nhzen sls linkv k chi fom foh We Br A_star ...
    B_star Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star ...
    De deltaY yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
    D_Matrix_Tm DD_Matrix_Tm tspan_star NTM rho R0 fung fung2 fungexp fungnlvis};

opts.POSDEF = true;
opts.SYM = true;

% Pre-allocating memory:
x = zeros(N,q,n+1);


t = linspace(0,tspan,n+1);
x(:,:,1) = xi;

xmin1(:,1) = min(x(:,:,1),[],2);
xmax1(:,1) = max(x(:,:,1),[],2);

%delta = 1.005; %inflation param
delta = beta;
%epsilon = 0.001;
l = peak_time - 1; % size of data assimilation window
Beta = (1/l).*ones(1,l); % MDA coefs (equal here)

% Special case: deleting points around first and second collapse
%{
num_del_pts = 1; % num points to delete around collapse (can be 0 for 
                 % only collapse point deleted)
if num_peaks > 2
    collapse_1_idx = peak_indices(2)-1;
else
    collapse_1_idx = peak_indices(1)-1;
end
Beta(collapse_1_idx-num_del_pts:collapse_1_idx+num_del_pts) = 0;
Beta(l-num_del_pts:l) = 0;
%}
%R = std^2;
timestep_time = zeros(1,n-l+1);
time_index = 1;
for j = 1
    %j % print index for debugging if it crashes
    %timestep_time(j) = toc;
    x10 = mean(x(:,:,j),2);
    A10 = squeeze(x(:,:,j)) - x10*ones(1,q);
    x1 = x10;
    TT = eye(q);
    
    dx1 = 1; %initialize greater than epsilon
    jj = 1;
    while norm(dx1) > epsilon & jj <= max_iter
        A1 = A10*TT;
        E1 = x1*ones(1,q) + A1;
        
        %E1(3,:) = max(E1(3,:),0.001); % Limit P to physical values
        
        TTinv = linsolve(TT,eye(size(TT)));
        
        t1 = t(time_index)/t0;
        t2 = t(time_index+1)/t0;
        
        parfor memb = 1:q
            %memb
            %tau_del{memb}
            [E2(:,memb,1),tau_del{memb}] = f(t1,t2,E1(:,memb),vars,tau_del{memb});
            %tau_del{memb}
        end
        y2(:,:,1) = h(E2(:,:,1));
        y2b(:,1) = mean(y2(:,:,1),2);
        HA2(:,:,1) = y2(:,:,1) - y2b(:,1)*ones(1,q);
        HA2(:,:,1) = HA2(:,:,1)*TTinv;
        dy2(:,1) = yth(:,time_index+1) - y2b(:,1);
        
        for kk = 2:l
            %kk
            t1 = t(time_index+kk-1)/t0;
            t2 = t(time_index+kk)/t0;
            for memb = 1:q
                [EE(:,memb),tau_del{memb}] = f(t1,t2,E2(:,memb,kk-1),vars,tau_del{memb});
            end
            E2(:,:,kk) = EE;
            y2(:,:,kk) = h(E2(:,:,kk));
            y2b(:,kk) = mean(y2(:,:,kk),2);
            HA2(:,:,kk) = y2(:,:,kk) - y2b(:,kk)*ones(1,q);
            HA2(:,:,kk) = HA2(:,:,kk)*TTinv;
            dy2(:,kk) = yth(:,time_index+kk) - y2b(:,kk);
            if disp_timesteps == 1
                disp(['iteration ',num2str(jj),', timestep ',num2str(kk),'/',num2str(l)])
            end
        end
        
        Rinv = linsolve(R,eye(size(R)));
        
        for kk = 1:l
            aux1(:,:,kk) = (HA2(:,:,kk)'*Beta(kk)*Rinv*HA2(:,:,kk))/(q-1);
            aux2(:,kk) = (HA2(:,:,kk)'*Beta(kk)*Rinv*dy2(:,kk))/(q-1);
        end
        
        
        GGinv = eye(q) + sum(aux1,3);
        GG = linsolve(GGinv,eye(size(GGinv)));
        b = linsolve(GGinv,sum(aux2,2));
        
        dx1 = A10*b + A10*GG*pinv(A10'*A10)*A10'*(x10-x1);
        %norm(dx1)
        if norm(dx1) <= epsilon | jj == max_iter
            x2 = mean(E2(:,:,end),2);
            A2 = E2(:,:,end) - x2*ones(1,q);
            E2a = x2*ones(1,q) + A2*delta;
        end
        
        x1 = x1 + dx1;
        TT = sqrtm(GG);
        
        x_est(:,jj) = x1;
        timestep_time(jj) = toc;
        %tracking_plots_4Dvar;
        disp(['ended iteration ',num2str(jj),' with norm(dx1) = ', ...
            num2str(norm(dx1)),' at ',num2str(timestep_time(jj)),' seconds'])
        jj = jj+1;
    end
    
    x(:,:,j) = E1;
%    x(:,:,j+1) = E2a;
    time_index = time_index + l;
    
end

c = clock;
%save(['IEnKS_MDA_lag',num2str(l),'_q',num2str(q),'test1_',datestr(now)],'-v7.3')
run_time = toc
