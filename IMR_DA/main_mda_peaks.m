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
%import_data
if strcmp(data_type,'sim') == 1
    import_data_sim
elseif strcmp(data_type,'sim_file') == 1
    load(data_filename)
elseif strcmp(data_type,'exp_Selda') == 1
    import_data_exp_Selda
else
    import_data_exp
end

% Create Initial Ensemble
create_ensemble_exp

tau_del = cell(q,1);
tau_del_f = tau_del;
%% Iterate

vars = {NT Pext_type Pext_Amp_Freq disptime Tgrad Tmgrad ...
    Cgrad comp t0 neoHook nhzen sls linkv k chi fom foh We Br A_star ...
    B_star Rv_star Ra_star L L_heat_star Km_star P_inf T_inf C_star ...
    De deltaY yk deltaYm xk yk2 Pv REq D_Matrix_T_C DD_Matrix_T_C ...
    D_Matrix_Tm DD_Matrix_Tm tspan_star NTM rho R0 fung fung2 fungexp fungnlvis};

opts.POSDEF = true;
opts.SYM = true;

% Pre-allocating memory:
x2 = zeros(N,q,n+1);


t = linspace(0,tspan,n+1);
x(:,:,1) = xi;

xmin1(:,1) = min(x(:,:,1),[],2);
xmax1(:,1) = max(x(:,:,1),[],2);

delta = beta;
Beta = (1/l).*ones(1,l); % MDA coefs (equal here)

n = peak_time - 1;
timestep_time = zeros(1,n-l+1);

for j = 1:(n-l+1)
    timestep_time(j) = toc;
    x10 = mean(x(:,:,j),2);
    A10 = squeeze(x(:,:,j)) - x10*ones(1,q);
    x1 = x10;
    TT = eye(q);
    
    dx = 1; % initialize greater than epsilon
    jj = 1;
    while norm(dx) > epsilon & jj <= max_iter
        A1 = A10*TT;
        E1 = x1*ones(1,q) + A1;
        
        TTinv = linsolve(TT,eye(size(TT)));
        
        t1 = t(j)/t0;
        t2 = t(j+1)/t0;
        
        tau_del = tau_del_f;
        parfor memb = 1:q
            [E2(:,memb,1),tau_del{memb}] = f(t1,t2,E1(:,memb),vars,tau_del{memb});
        end
        y2(:,:,1) = h(E2(:,:,1));
        y2b(:,1) = mean(y2(:,:,1),2);
        HA2(:,:,1) = y2(:,:,1) - y2b(:,1)*ones(1,q);
        HA2(:,:,1) = HA2(:,:,1)*TTinv;
        dy2(:,1) = yth(:,j+1) - y2b(:,1);
        
        for kk = 2:l
            t1 = t(j+kk-1)/t0;
            t2 = t(j+kk)/t0;
            parfor memb = 1:q
                [EE(:,memb),tau_del{memb}] = f(t1,t2,E2(:,memb,kk-1),vars,tau_del{memb});
            end
            E2(:,:,kk) = EE;
            y2(:,:,kk) = h(E2(:,:,kk));
            y2b(:,kk) = mean(y2(:,:,kk),2);
            HA2(:,:,kk) = y2(:,:,kk) - y2b(:,kk)*ones(1,q);
            HA2(:,:,kk) = HA2(:,:,kk)*TTinv;
            dy2(:,kk) = yth(:,j+kk) - y2b(:,kk);
        end
        
        Rinv = linsolve(R,eye(size(R)));
        
        for kk = 1:l
            aux1(:,:,kk) = (HA2(:,:,kk)'*Beta(kk)*Rinv*HA2(:,:,kk))/(q-1);
            aux2(:,kk) = (HA2(:,:,kk)'*Beta(kk)*Rinv*dy2(:,kk))/(q-1);
        end
        
        
        GGinv = eye(q) + sum(aux1,3);
        GG = linsolve(GGinv,eye(size(GGinv)));
        b = linsolve(GGinv,sum(aux2,2));
        
        dx = A10*b + A10*GG*pinv(A10'*A10)*A10'*(x10-x1);
        
        if norm(dx) <= epsilon | jj == max_iter
            x2 = mean(E2(:,:,1),2);
            A2 = E2(:,:,1) - x2*ones(1,q);
            E2a = x2*ones(1,q) + A2*delta;
            disp(['timestep ',num2str(j),'/',num2str(n-l+1), ...
                ': norm(dx1) = ',num2str(norm(dx)),', ',num2str(jj), ...
                ' iterations'])
            
            tau_del_f = tau_del;
        end
        
        x1 = x1 + dx;
        TT = sqrtm(GG);
        
        jj = jj+1;
    end
   
    x(:,:,j) = E1;
    x(:,:,j+1) = E2a;
    
end

c = clock;
run_time = toc