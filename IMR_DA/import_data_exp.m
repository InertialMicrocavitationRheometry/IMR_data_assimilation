% This file exports data in vector yth, it must be re-written for each data
% set to ensure the data is formatted correctly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Data format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_data = load([data_filepath,data_filename]); % load file 
Rnew = exp_data.M(:,2); % extract R values from loaded file (normalized, unscaled)
t = exp_data.M(:,1); % extract t values from loaded file (normalized, unscaled)

if length(Rnew(1,:)) == 1
    Rnew = Rnew';
end
if length(t(1,:)) == 1
    t = t';
end

Rexp = Rnew;%*1e-6;
[R0,max_index] = max(Rexp);
%yth = Rnew;
yth = Rexp(max_index:end)./R0;
t = t*R0/sqrt(101325/1060);
t = t(max_index:end);

% plot(t,yth,'k')
plot(1:length(yth),yth,'rx','linewidth',2)

kk = 1;
for jj = 1:length(yth)
    if isnan(yth(jj))
    else
        yth(kk) = yth(jj);
        t(kk) = t(jj);
        kk = kk + 1;
    end
end
yth = yth(1:kk-1);
t = t(1:kk-1);

tspan = t(end)-t(1)
n = length(t)-1;

if exist('l') == 0
    l = n;
end
timesteps = n+1;
% Find peak_time
peak_indices = find(islocalmin(yth));
peak_time = peak_indices(num_peaks);