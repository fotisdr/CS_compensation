% Script to extract, analyse, and plot EFR responses from recorded bdf files
% The names of the bdf files need to be defined in the condition array
% 
% Fotis Drakopoulos, February 2021, UGent
%

clc
clear
% close all

data_dir = ''; % folder with the bdf files to analyse
save_dir = ['Analysis/']; % folder to save the analysis results

% names of the bdf files (without the extension)
condition=[{'SAM_4kHz_120Hz'}, {'David_220Hz_hp'}];

%% Define the parameters for the analysis

opts.FS=16384; % sampling frequency
% opts.f0=120; % fundamental frequency of the stimulus - for plotting purposes
opts.f0=[120,220]; % can be also defined separately for every condition

% duration of stimulus in secs
opts.stim_dur = [0.5,0.5];
% secs to omit from the onset of the stimulus or to add before the trigger 
opts.onset_dur = [0.1,0.00005];

% opts.threshold=25; % epoch rejection threshold (p-p) - if not set a dynamic threshold is used

opts.channels=48; % channels to analyse
% 48 is the Cz, [11,12,19,32,46:49,56] for multi-channel around Cz
opts.CR1=65; % reference channel 1
opts.CR2=66; % reference channel 2

% Filter the EEG responses - for more information refer to the EEG_Analysis script
opts.filtering = 1; % 0 applies no filtering, 1 applied bandpass filtering

% Apply bootstrap to the recorded EEG signals - for more information refer to the EEG_Analysis script
opts.bootstrap = 1;         % 0 or 1 to apply or not bootstrapping to the signal
opts.save_bootstrapped = 1; % save bootstrapped signal
opts.plot = 1;              % plot results

%% Stage 1: Extract bdf to mat
addpath('functions')
channels = opts.channels;
chunks = 100;

% check if files are already extracted
for j=1:length(channels)
    for i=1:length(condition)
        if ~exist([save_dir filesep condition{i} 'Ch' num2str(channels(j)) '.mat'],'file')
            extract_flag(i,j)=1;
        else
            extract_flag(i,j)=0;
        end
    end
end

% extract bdf to mat
for i=1:length(condition)
    if all(extract_flag(i,:))
        disp('Stage 1: Extract bdf to mat');
        extract_bdf_data([data_dir filesep condition{i} '.bdf'],[],[],chunks,'off')
    end
end

%% Stage 2: Load triggers
for i=1:length(condition)
    if any(extract_flag(i,:))
        if ~exist([save_dir filesep],'dir')
            mkdir([save_dir filesep])
        end
        disp('Stage 2: Load triggers');
        TriggerLoader([data_dir filesep condition{i}],[save_dir filesep],condition{i})
    end
end

%% Stage 3: Load channel data
if any(extract_flag(:))
    disp('Stage 3: Load channel data');
    all_channels_loading(condition,opts.channels,opts.CR1,opts.CR2,data_dir,save_dir)
end

%% Stage 4: EEG analysis
opts.analysis_dur = opts.stim_dur - opts.onset_dur; % total_stimulus_duration - onset_duration
output = EEG_Analysis(condition,opts,save_dir);
