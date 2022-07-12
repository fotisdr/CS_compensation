function [] = run_experiment(stim,pars,cs)
% This function generates and plays a stimulus in MATLAB using Playrec
% 
% A stim array must be given as input
% The second argument (pars) is a struct array that contains the recording 
% setup parameters
% A number can be given as the third argument (cs), corresponding to the 
% respective CS profile, to process the stimulus before playback
% 
% Adjusted by Fotis Drakopoulos, December 2019, UGent
%

%% Set the needed paths for the experiment
if nargin < 2
    pars.device_play = 0; % Device to play
    pars.device_rec = 0; % Device to record
    pars.play_chan = [1,9]; % Channels to use for playing
    pars.rec_chan = [1,9]; % Channels to use recording

    pars.dur = 0.5; % duration of epoch in sec
    pars.avs = 1000; % repetions
    pars.jit = 0.1; % percentage of the silence gap
    pars.silence = 100e-3; %silence gap in sec
    pars.fs = 48000; % stimulus sampling rate
    pars.loops = 10; % number of loops to devide the stimuli in
end

% uncomment for the two stimuli provided
%%%%%%%%%%%% SAM stimulus (SAM_4kHz_120Hz.mat)
% pars.fs = 48000
% pars.loops = 10
% pars.reps = 1000
% pars.win_len = 0; % for processing the SAM tones
%%%%%%%%%%%% Speech-Matrix stimulus (David_220Hz_hp.mat)
% pars.fs = 44100
% pars.loops = 6
% pars.reps = 1200
% % parameters for applying processing on speech
% pars.win_len = 25; % length of the running RMS window 
% pars.height = 0; % minimum height of envelope peaks to be considered
% pars.prominence = 1.15 % minimum difference from adjacent peaks
% pars.distance = 1; % minimum distance between peaks

if nargin < 3
    cs = 1333;
else
    addpath('../algorithms');
end

tic
% Make sure we can see external m files, even if this file is not in the
% current working directory but on the path somewhere
MAINPATH = mfilename('fullpath'); % get the full path of this function (+ its name)
MAINPATH = MAINPATH(1:end-length(mfilename)); % substract the name from the path

addpath(MAINPATH);
addpath([MAINPATH,'functions']); % Add path were all other needed files are located + playrec mex32

stim=stim(1:round(pars.dur*pars.fs));

% use three ones instead of one to make audio and eeg samplingrate work
% together
TrigEpoch=[1,1,1, zeros(1,ceil(pars.dur*pars.fs)+ceil(pars.fs*(pars.silence+pars.jit*pars.silence)))];

seed = sum(uint8(stim)); % seed for the rng; different for ev                                            ery condition --> stimname expressed as numbers and summed up
rng(seed,'twister') %sets the random no generator to a constant

JitN=round(pars.fs*2*pars.jit*pars.silence*rand(pars.avs,1)); %creates all jitter instances

if rem(pars.avs,pars.loops) ~= 0
    error('Number of repetions is not a multiple of the number of loops')
end
if rem(pars.avs,2) ~= 0
    error('Number of both polarities is unequal; set pars.avs to multiple of two')
end

%% Initialise playrec and play sound

if playrec('isInitialised') == 1   % tests if playrec is initialised allready - perhabs earlier using did not reset it
    playrec('reset')    ;          % in state of initialisation playrec can't be initialised, therefore reset it before using
end

playrec('init',pars.fs,pars.device_play,pars.device_rec);  % initialisation of playrec


for z = 1:pars.loops % for the number of loops
    
    Stim=[]; % initiate stimulus
    Trig=[]; % initiate trigger stream
    
    fprintf(':: Starting loop %d with %d epochs\n\n', z,pars.avs/pars.loops)
    
    for k=1:pars.avs/pars.loops % for the amount of epochs per loop
        
        if rem(k,2)==1 % for odd epochs
            modnoise = stim;
        else
            modnoise=-modnoise;
        end
        
        Epoch = [modnoise zeros(1,round(pars.fs*(pars.silence+pars.jit*pars.silence)))]; 
        Stim = [Stim Epoch(1:numel(Epoch)-JitN(k))]; 
        Trig = [Trig TrigEpoch(1:numel(Epoch)-JitN(k))];
        
    end
    
    if cs ~= 1333 % process the generated stimulus
        if win_len % RMS-based envelope estimation
            Stim = g_70dB(Stim,cs,win_len,height,prominence,distance);
        else % Hilbert envelope estimation
            Stim = g_70dB(Stim,cs);
        end
    end

    prepostPause = [0.1 0.1]; % Amount of silence at the begiining and and of a loop in seconds
    
    Ch1 = [zeros(1,round(pars.fs*prepostPause(1)))  Stim  zeros(1,round(pars.fs*prepostPause(2)))];
    Trigger = [zeros(1,round(pars.fs*prepostPause(1)))  Trig  zeros(1,round(pars.fs*prepostPause(2)))];
    trigcode = bin2dec('10000000 0000 0000 0000 0000')/(1*16777215); % trigger code
    Trigger = Trigger*trigcode; % turn trigger ones into trigger code
    
%     Ch1 = Ch1*10^((pars.level_total-pars.lcal_total )/20); %here the amplitude is set
%     fprintf(':: Calculated SPL (re: 20uPa) of loop: %0.2f\n', 20*log10(rms(Ch1)/20e-6))
    
  rms(Ch1)
    if max(abs(Ch1)) >= 1
        warning('Outputlevel exeeds 1; this might lead to clipping of the stimulus')
    end
    
    PlayBuffer = [Ch1;Trigger]'; % left, right
    
    disp(':: recording started - please wait until finished')
    
    playrec_buf(PlayBuffer,pars.device_play, pars.device_rec, pars.play_chan, pars.rec_chan, pars.fs); % Main function that playes sound
    
end

playrec('reset') ;  % closing the initialisation
clear('playrec') ;  % clearing all playrec-memories


disp(':: Condition has finished')

toc
end



