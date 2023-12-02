% This file exports data in vector yth, it must be re-written for each data
% set to ensure the data is formatted correctly

%dataset = 8; %between 1-20 for water case
%addpath ./IMR-vanilla/functions
%addpath ./IMR-vanilla/src
% load('/scratch/jspratt/EnKS/exp_data/11kPa_PA/RofTdata.mat')
%% old data format
%{
load([data_filepath,data_filename]);
Rexp = Rnew(dataset,:)*1e-6;
[R0,max_index] = max(Rexp);

yth = Rexp(max_index:end)./R0;
t = t(max_index:end);

tspan = t(end)-t(1);
n = length(t)-1;

if exist('l') == 0
    l = n;
end
timesteps = n+1;

% Find peak_time
peak_indices = find(islocalmin(yth));
peak_time = peak_indices(num_peaks);
%}

%% new data format
exp_data = load([data_filepath,data_filename]);
%Rnew = exp_data.R;
Rnew = exp_data.PolyAcry_12_48_08(:,2)';
%t = exp_data.t;
t = exp_data.PolyAcry_12_48_08(:,1)';

% fixing shape if needed
if length(Rnew(1,:)) == 1
    Rnew = Rnew';
end
if length(t(1,:)) == 1
    t = t';
end

Rexp = Rnew*1e-6;
[R0,max_index] = max(Rexp);

yth = Rexp(max_index:end)./R0;
t = t(max_index:end);

% delet NaNs
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
%

tspan = t(end)-t(1);
n = length(t)-1;

if exist('l') == 0
    l = n;
end
timesteps = n+1;

% Find peak_time
peak_indices = find(islocalmin(yth));
peak_time = peak_indices(num_peaks);

% Exceptionally, for test run without collapse points
%{
if method == 'EnKS'
    
    if num_peaks > 2 % checking for indexing error present in some runs
        collapse_1_idx = peak_indices(2);
    else
        collapse_1_idx = peak_indices(1);
    end
    
    yth = [yth(1:collapse_1_idx-2),yth(collapse_1_idx+2:peak_time-2), ...
        yth(peak_time+2:end)];
    t = [t(1:collapse_1_idx-2),t(collapse_1_idx+2:peak_time-2), ...
        t(peak_time+2:end)];
    
    tspan = t(end) - t(1);
    n = length(t) - 1;
    timesteps = n+1;
end
%}