% This is the plotting file for the data assimilation paper, it is split up
% into different sections for different data types, methods, and type of
% plots. It is meant to be run locally
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  add paths, define variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextinterpreter','latex');
% For dimensional parameter calculation:
rho = 1060; % (Kg/m^3) Liquid Density taken from Estrada et. al
Uc = sqrt(P_inf/rho);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% extract needed data from ensemble:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j_range = 1:kk+1;
jj_range = 0:jj-1;

Ca_est = exp(x_est(2*NT+NTM+7,1:jj-1));
Re_est = exp(x_est(2*NT+NTM+8,1:jj-1));
De_est = x_est(2*NT+NTM+9,1:jj-1);
alpha_est = exp(x_est(2*NT+NTM+10,1:jj-1));
lambda_nu_est = x_est(2*NT+NTM+11,1:jj-1);

Ca_est = [Ca,Ca_est];
Re_est = [Re,Re_est];
De_est = [De,De_est];
alpha_est = [alpha,alpha_est];
lambda_nu_est = [lambda_nu,lambda_nu_est];

%CaError = 100*(abs(Ca_est-Ca_true)./Ca_true);
%ReError = 100*(abs(Re_est-Re_true)./Re_true);

G_est = P_inf./Ca_est;
%GError = 100*(abs(G_est-G_true)./G_true);

mu_est = (P_inf*R0)./(Re_est.*Uc);
%muError = 100*(abs(mu_est-mu_true)./mu_true);

%runlims = timesteps:timesteps:iterations*timesteps; % specify limits of runs for tables and plots
runlims = j_range(end);
xf = [x1,squeeze(mean(E2,2))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Ensemble statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = cat(3,xi,E2);
x_min = squeeze(min(E,[],2));
x_max = squeeze(max(E,[],2));
clear std; %interferes with matlab std function
sigma1 = squeeze(std(E,0,2));
x_minsigma = xf-sigma1;
x_plussigma = xf+sigma1;

% stats for G,mu (through Ca,Re)
X = cat(3,xi,x(:,:,1:jj-1));

X_min = squeeze(min(X,[],2));
X_max = squeeze(max(X,[],2));
sigma2 = squeeze(std(X,0,2));
X_minsigma = squeeze(mean(X,2))-sigma2;
X_plussigma = squeeze(mean(X,2))+sigma2;

Ca_min = exp(X_min(2*NT+NTM+7,:));
Ca_max = exp(X_max(2*NT+NTM+7,:));
Ca_minsigma = exp(X_minsigma(2*NT+NTM+7,:));
Ca_plussigma = exp(X_plussigma(2*NT+NTM+7,:));

Re_min = exp(X_min(2*NT+NTM+8,:));
Re_max = exp(X_max(2*NT+NTM+8,:));
Re_minsigma = exp(X_minsigma(2*NT+NTM+8,:));
Re_plussigma = exp(X_plussigma(2*NT+NTM+8,:));

De_min = exp(X_min(2*NT+NTM+9,:));
De_max = exp(X_max(2*NT+NTM+9,:));
De_minsigma = exp(X_minsigma(2*NT+NTM+9,:));
De_plussigma = exp(X_plussigma(2*NT+NTM+9,:));

alpha_min = exp(X_min(2*NT+NTM+10,:));
alpha_max = exp(X_max(2*NT+NTM+10,:));
alpha_minsigma = exp(X_minsigma(2*NT+NTM+10,:));
alpha_plussigma = exp(X_plussigma(2*NT+NTM+10,:));

lambda_nu_min = exp(X_min(2*NT+NTM+11,:));
lambda_nu_max = exp(X_max(2*NT+NTM+11,:));
lambda_nu_minsigma = exp(X_minsigma(2*NT+NTM+11,:));
lambda_nu_plussigma = exp(X_plussigma(2*NT+NTM+11,:));

G_min = P_inf./Ca_min;
G_max = P_inf./Ca_max;
G_minsigma = P_inf./Ca_minsigma;
G_plussigma = P_inf./Ca_plussigma;

mu_min = (P_inf*R0)./(Re_min.*Uc);
mu_max = (P_inf*R0)./(Re_max.*Uc);
mu_minsigma = (P_inf*R0)./(Re_minsigma.*Uc);
mu_plussigma = (P_inf*R0)./(Re_plussigma.*Uc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Other metrics to plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_err = abs(xf(1,:)-yth(1:time_index));
avg_R_err = mean(R_err);
R_percent_err = 100*(R_err./yth(1:time_index));
avg_R_percent_err = mean(R_percent_err);
RMSE = sqrt(mean(R_err.^2));
NRMSE = RMSE/(mean(yth(1:time_index)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R evolution with final initial conditions
figure(1)
clf
hold on
%plot(1:length(xth(1,:)),xth(1,:),'k-','linewidth',2)
plot(1:length(yth),yth,'rx','linewidth',2)
plot(1:length(xf(1,:)),xf(1,:),'bo-')
fill([j_range,fliplr(j_range)],[x_min(1,j_range),fliplr(x_max(1,j_range))], ...
    'b','FaceAlpha',0.2,'EdgeColor','none')
fill([j_range,fliplr(j_range)],[x_minsigma(1,j_range),fliplr(x_plussigma(1,j_range))], ...
    'b','FaceAlpha',0.1,'EdgeColor','none')
xlim([0 peak_time])
legend('measurement','4Dvar estimate','location','northwest')
xlabel('time step')
ylabel('$\frac{R}{R_{max}}$')
grid on
set(gca,'fontsize',20)


% % U, P, S evolution with final initial conditions
% 
% figure(2)
% clf
% subplot(3,1,1)
% hold on
% %plot(1:length(xth(2,:)),xth(2,:),'k-','linewidth',2)
% plot(1:length(xf(2,:)),xf(2,:),'bo-')
% fill([j_range,fliplr(j_range)],[x_min(2,j_range),fliplr(x_max(2,j_range))], ...
%     'b','FaceAlpha',0.2,'EdgeColor','none')
% fill([j_range,fliplr(j_range)],[x_minsigma(2,j_range),fliplr(x_plussigma(2,j_range))] ...
%     ,'b','FaceAlpha',0.1,'EdgeColor','none')
% %xlim([0 50])
% legend('4Dvar estimate','location','northwest')
% ylabel('U')
% grid on
% set(gca,'fontsize',20)
% 
% subplot(3,1,2)
% hold on
% %plot(1:length(xth(3,:)),xth(3,:),'k-','linewidth',2)
% plot(1:length(xf(4,:)),exp(xf(3,:)),'bo-')
% fill([j_range,fliplr(j_range)],[exp(x_min(3,j_range)),fliplr(exp(x_max(1,j_range)))], ...
%     'b','FaceAlpha',0.2,'EdgeColor','none')
% fill([j_range,fliplr(j_range)],[exp(x_minsigma(3,j_range)),fliplr(exp(x_plussigma(1,j_range)))] ...
%     ,'b','FaceAlpha',0.1,'EdgeColor','none')
% %xlim([0 50])
% ylabel('$\frac{p_b}{p_{\infty}}$')
% grid on
% set(gca,'fontsize',20)
% 
% subplot(3,1,3)
% hold on
% %plot(1:length(xth(4,:)),xth(4,:),'k-','linewidth',2)
% plot(1:length(xf(4,:)),xf(4,:),'bo-')
% fill([j_range,fliplr(j_range)],[x_min(4,j_range),fliplr(x_max(4,j_range))], ...
%     'b','FaceAlpha',0.2,'EdgeColor','none')
% fill([j_range,fliplr(j_range)],[x_minsigma(4,j_range),fliplr(x_plussigma(4,j_range))] ...
%     ,'b','FaceAlpha',0.1,'EdgeColor','none')
% %xlim([0 50])
% xlabel('time step')
% ylabel('$\frac{S}{p_{\infty}}$')
% grid on
% set(gca,'fontsize',20)
% 
% % log(Ca), log(Re)
% figure(8)
% clf
% 
% subplot(2,1,1)
% hold on
% plot(jj_range,log(Ca_est),'ko-')
% legend('Estimate','location','northeast')
% xlabel('iteration')
% ylabel('$log(Ca)$')
% grid on
% set(gca,'fontsize',20)
% 
% subplot(2,1,2)
% hold on
% plot(jj_range,log(Re_est),'ko-')
% xlabel('iteration')
% ylabel('$log(Re)$')
% grid on
% set(gca,'fontsize',20)
% 
% % Compare truth and estimate for G,mu:
% figure(89)
% clf
% subplot(2,1,1)
% hold on
% %plot(jj_range,G_true*ones(1,jj),'k-','linewidth',2)
% plot(jj_range,G_est,'ko-')
% %fill([jj_range,fliplr(jj_range)],[G_min,fliplr(G_max)],'b','FaceAlpha',0.2,'EdgeColor','none')
% %fill([jj_range,fliplr(jj_range)],[G_minsigma,fliplr(G_plussigma)],'b','FaceAlpha',0.1,'EdgeColor','none')
% legend('Estimate','location','northwest')
% ylabel('$G \ [Pa]$')
% grid on
% set(gca,'fontsize',20)
% 
% subplot(2,1,2)
% hold on
% %plot(jj_range,mu_true*ones(1,jj),'k-','linewidth',2)
% plot(jj_range,mu_est,'ko-')
% %fill([jj_range,fliplr(jj_range)],[mu_min,fliplr(mu_max)],'b','FaceAlpha',0.2,'EdgeColor','none')
% %fill([jj_range,fliplr(jj_range)],[mu_minsigma,fliplr(mu_plussigma)],'b','FaceAlpha',0.1,'EdgeColor','none')
% xlabel('iteration')
% ylabel('$\mu \ [Pa \cdot s]$')
% grid on
% set(gca,'fontsize',20)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% Other parameters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(90)
% clf
% subplot(2,1,1)
% hold on
% plot(jj_range,alpha_est,'ko-')
% legend('Estimate','location','northeast')
% xlabel('iteration')
% ylabel('$\alpha)$')
% grid on
% set(gca,'fontsize',20)
% 
% subplot(2,1,2)
% hold on
% plot(jj_range,lambda_nu_est,'ko-')
% xlabel('iteration')
% ylabel('$\lambda_{\nu}$')
% grid on
% set(gca,'fontsize',20)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Paper plots June 29
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Radius with normalized time
% figure(111)
% clf
% hold on
% %plot(1:length(xth(1,:)),xth(1,:),'k-','linewidth',2)
% plot(t(j_range)./t0,yth(j_range),'rx','linewidth',2)
% plot(t(j_range)./t0,xf(1,j_range),'bo-')
% fill([t(j_range)./t0,fliplr(t(j_range)./t0)],[x_min(1,j_range),fliplr(x_max(1,j_range))], ...
%     'b','FaceAlpha',0.2,'EdgeColor','none')
% fill([t(j_range)./t0,fliplr(t(j_range)./t0)],[x_minsigma(1,j_range),fliplr(x_plussigma(1,j_range))], ...
%     'b','FaceAlpha',0.1,'EdgeColor','none')
% xlim([0 t(length(j_range))./t0])
% legend('measurement','4Dvar estimate','ensemble','location','northwest')
% xlabel('$t^*$')
% ylabel('$\frac{R}{R_{max}}$')
% grid on
% set(gca,'fontsize',20)
% %}