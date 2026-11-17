%%% Optimisation of PUSH pulses for 2D or 3D imaging

clearvars; close all; clc; 

% Refer to https://github.com/mriphysics/push for dependent PUSH scripts
addpath(genpath('/insert_your_path_to_scripts/PUSH/'))
% Need to contact Siemens for proprietary code
addpath('/insert_your_path_to_siemems_script_dependencies/siemens')

if isempty(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = 10;
    parpool(c, c.NumWorkers);
end
tic 

%% paths to fieldmaps for UP
% use make_AdjMask.m to make mask.mat

% Adjustment data from multiple participants
file_paths = {'/insert_your_path_to_corresponding_fieldmaps1/AdjDataUser',
              '/insert_your_path_to_corresponding_fieldmaps2/AdjDataUser',
              '.....'
              };
% Brain masks for respective participants
mask_paths = {'/insert_your_path_to_corresponding_mask1/mask.mat',
              '/insert_your_path_to_corresponding_mask2/mask.mat',
              '.....'
              };


% combine Adjustment data

tx = []; % list of B1 maps for each subject
mask = []; % list of masks for each subject
for i = 1:length(file_paths)
    file_i = load(strcat(file_paths{i}, '.mat'));
    dims = size(file_i.Adj.B0);
    Nr = prod(dims); slices = 1:dims(3);
    Nch = file_i.Adj.coils;
    S_i = file_i.Adj.S;  % coil sensitivity map of subject i

    load(mask_paths{i});%file_i.Adj.W;
    mask_i = mask_final(:,:,slices); % brain mask for subject i

    tx_i = permute(S_i,[1 3 2]);
    txmaps_i = reshape(tx_i,[dims Nch]);
    tx_i = reshape(txmaps_i,[Nr Nch]); 

    tx = [tx; tx_i];
    mask = [logical(mask);mask_i(:)];
end

%% Load and pre-process tx maps
 
% load('AdjDataUser.mat')
Sar = load('SarDataUser.mat');
Sys = load('SysDataUser.mat');     
Sys.CPDvalid = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set optimisation options and define (SAR & hardware) constraints

opt_options.w0                = [];
opt_options.Niter             = 1e3;
opt_options.quality           = 4;
opt_options.verbose           = 'off'; %{'off', 'final', 'iter}
opt_options.multistart        = true;
opt_options.multistart_trials = 1e1;
opt_options.weights           = [];
opt_options.rng_seed          = 'default';

%TR  = 22e-3; %[s]
Nburst = 10; % for MT boost, No. of bursts
TR_burst = 2000e-3 / Nburst;
tau = 4e-3;  %[s]
dt  = 10e-6; %[s] dwell time
b1waveform = gen_B1_pulse(1, tau, dt, 'Gaussian'); b1waveform = b1waveform ./ max(abs(b1waveform));
avg2peak = sqrt( TR_burst / (sum(abs(b1waveform).^2) * dt)); %scaling from average B1rms to peak B1rms

target_b1rms_pulse = 4.5; % uT

%%% add a global-SAR equivalent VOP (type 99, pulse energy)
if isfield(Sys,'Y3Matrix') && ~isempty(Sys.Y3Matrix) && ~isempty(Sar)
    ind=find(Sar.ZZtype==5); 
    if length(ind)~=1
        cpl_vop = zeros(size(Sys.Y3Matrix));
    else
        cpl_vop = Sar.ZZ(:,:,ind);
    end
    Sar.ZZ(:,:,end+1) = Sys.Y3Matrix - cpl_vop;
    Sar.ZZtype(end+1) = 99;
end

opt_options.power_factor_s = sum(abs(b1waveform).^2)*dt; %for SAR calculations
opt_options.VOP_type = [1; 2; 6; 7; 8; 99];

% this line chooses which VOPs are relevant / used 
% VOP: 6, 8, and "99"
% 6 = Local SAR
[~,~,SEDtype] = calc_SED(Sys,Sar,ones(Nch,1),Nch,1e-3,opt_options.VOP_type,'100');

SED_limit = zeros(size(SEDtype));
for h = 0:max(SEDtype)
    if h~=99
        SED_limit(SEDtype==h) = Sar.VOPTypeLimit(Sar.VOPTypeLimitID==h);% * (TR_ms/1000);
    else
        SED_limit(SEDtype==h) = Sar.GlobalSARLimit / Sar.AF;%* (TR_ms/1000) / Sar.AF;
    end
end
% Account for sequence timing
SED_limit = SED_limit * TR_burst;

%if not in first level mode, apply it
if (Sar.VOPTypeLimitID(2)==6 && Sar.VOPTypeLimit(2)==10)
    SED_limit(SEDtype==6) = 2 * SED_limit(SEDtype==6); 
end
opt_options.SED_limit = SED_limit;


%% Optimisation

filename = ['./bin/PUSH_opt_UP_20subj_Tr_', num2str(TR_burst*1e3), 'ms_gauss_', num2str(tau*1e3),'ms.mat']; 
    
subpulses = 2;                    
Nsp = numel(subpulses);

% target_avgB1rms = 0.8:0.1: 1.1; %0.8:0.2:2.0;    
target_avgB1rms = target_b1rms_pulse *sqrt(subpulses*tau/TR_burst); 
Ntg = numel(target_avgB1rms); 


all_wopt = cell(Nsp,Ntg); 
all_fval = zeros(opt_options.multistart_trials,Nsp,Ntg); 
all_runtime = zeros(Nsp,Ntg);

CPmode_peakB1rms = abs(tx * Sys.txScaleFactor); 
CPmode_wopt = cell(Ntg,1);

tx_mask = tx(mask(:),:);

rowsToKeep = ~all(tx_mask == 0, 2);

% Use this logical index to keep only the desired rows
tx_mask = tx_mask(rowsToKeep, :);

%par
for tt=1:Ntg

    % create copy of opt_options to allow multi-threading 
    % /!\ keep opt_options unchanged inside parfor /!\
    cp_opt_options = opt_options;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           CP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimal voltage for CP mode (minimises the cost function):
    cp_opt_options.base_voltage = (target_avgB1rms(tt) * avg2peak) * sum(CPmode_peakB1rms(mask)) / sum(CPmode_peakB1rms(mask).^2);
    CPmode_wopt{tt} = cp_opt_options.base_voltage * Sys.txScaleFactor;
    % % check if optimal CP voltage is within SED limits
    F_sed = calc_SED(Sys,Sar,CPmode_wopt{tt},Nch,opt_options.power_factor_s,opt_options.VOP_type,'100');
    relative_F_sed = F_sed ./ opt_options.SED_limit; 
    % if above SED limits, downweight CP voltage until it's feasible
    while max(relative_F_sed)>1
        CPmode_wopt{tt} = 0.9999 * CPmode_wopt{tt};
        F_sed = calc_SED(Sys,Sar,CPmode_wopt{tt},Nch,opt_options.power_factor_s,opt_options.VOP_type,'100');
        relative_F_sed = F_sed ./ opt_options.SED_limit; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   optimise PUSH pulses
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    target_peakB1rms = target_avgB1rms(tt) * avg2peak * ones(size( tx(mask(:),:),1),1);
    
    tic
    
    [all_wopt{tt}, all_fval(:,tt)] = optimise_PUSH_B1rms(tx_mask, target_peakB1rms(rowsToKeep), subpulses, Sys, Sar, cp_opt_options);
    all_runtime(tt) = toc;                                           
    
    fprintf(1,'Finished optimisations using B1rms target = %.2f uT.\n',target_avgB1rms(tt))
end
save(filename)
toc
