function output=EEG_Analysis(condition,opts,savedir)
% Script to analyse EEG data 
% The condition array defines the recorded conditions
% The opts struct includes all the analysis parameters
% The results are saved in savedir
%
% Originally created by Sarineh Keshishzadeh, UGent
% Adapted by Fotis Drakopoulos, UGent, April 2020
% 

if nargin < 3
    savedir = '';
else
    savedir = [savedir filesep];
end

FS = opts.FS; % sampling frequency
f0 = opts.f0; % fundamental frequency of the stimulus
% stim_dur = opts.stim_dur; % duration of stimulus in secs
onset_dur = opts.onset_dur; % secs to omit from the onset of the stimulus
analysis_dur = opts.analysis_dur; 

if exist('opts.threshold','var')
    Threshold=opts.threshold; % epoch rejection threshold (p-p)
end

channels = opts.channels; % channels to analyse
CR1 = opts.CR1; % reference channel 1
CR2 = opts.CR2; % reference channel 2

%% analyse each condition
for condi=1:length(condition)
    name=[savedir cell2mat(condition(condi))];
        
    %% Triggers
    load([savedir 'Trigs',cell2mat(condition(condi)),'.mat']); % load the triggers
    % unique(Trigs); % show values of triggers
    
    T=find(diff(Trigs)==1 | diff(Trigs)==2)+1; % Detect triggers of 1 or 2
    disp(['Found triggers: ' num2str(length(T))])

    %% Load data
    chani=0;
    for channo = channels
        chani = chani + 1;
        load([name,'Ch',num2str(channo),'.mat'])
        disp(['Loaded Channel ' num2str(channo) ' data'])
        
        %% EEG filtering
        if (length(opts.filtering)==length(condition) & opts.filtering(condi)==1) | ...
                (length(opts.filtering)==1 & opts.filtering==1)
            if opts.f0(condi) == 120
                fclo=600; %cut-off frequency (upper limit)
            else
                fclo=2500; %cut-off frequency (upper limit)
            end
            fchi=60; %cut-off frequency (lower limit)

            b = fir_bp(fchi,fclo,FS,800); % digital bandpass filtering
            Data=filtfilt(b,1,Ch);
        else
            Data=Ch;
        end

        %% Epoching happens here
        if length(onset_dur)==length(condition) & length(onset_dur)==length(condition)
            starti=round(onset_dur(condi)*FS);
            endi=round(analysis_dur(condi)*FS);
        else
            starti=round(onset_dur*FS);
            endi=round(analysis_dur*FS);
        end
        
        ER=zeros(numel(T),endi+1);
        for n=1:numel(T)
            ER(n,:)=Data(T(n)+starti:T(n)+starti+endi);
        end
        
        time_vec_rec = (starti:1:(size(ER,2)+starti-1))/FS;
        freq_vec_rec = 0:FS/(size(ER,2)):FS-FS/size(ER,2);
        
        %% Base line correction
        ER=ER-mean(ER,2);
        
        vmax = max(ER,[],2); % max voltage per epoch
        vmin = min(ER,[],2); % min voltage per epoch
        vpp = abs(vmax-vmin); % positive difference from most negative to postive voltage (p-

        if ~exist('Threshold','var')
            % compute a dynamic threshold based on the moving median
            [~,~,Thresh,~] = isoutlier(vpp,'movmedian',200,'ThresholdFactor',3);
            % Thresh=prctile(vpp,95);
        else
            Thresh=Threshold;
        end
        
        % visualise epoch rejection threshold
        if opts.plot
            figure
            plot(vpp), hold on
            plot(1:size(vpp), repmat(Thresh, 1, size(ER,1)), 'r--');
            legend('Data','Threshold'), xlabel('Samples');
            drawnow()
        end

        good = find(vpp < Thresh); % good epochs
        
        ERp = ER(good(mod(good,2)==1),:);  % positive polarity epochs
        ERn = ER(good(mod(good,2)==0),:);  % negative polarity epochs
        
        ER = ER(good,:); % reject bad epochs
        T = T(good);
        
        fprintf('\nEpochs preserved: %d\n',size(ER,1))

        % PP = ERp'; % positive polarity epochs
        % PN = ERn'; % negative polarity epochs
        
        %% EFR FFT
        NFFT = 2^nextpow2(size(ER,2));
        fprintf('FFT length: %d samples\n',NFFT)
        
        % zero-pad the time-domain signals
        ER = [ER zeros(size(ER,1),NFFT-size(ER,2))];
        ERp = [ERp zeros(size(ERp,1),NFFT-size(ERp,2))];
        ERn = [ERn zeros(size(ERn,1),NFFT-size(ERn,2))];
        
        freq=FS/2*linspace(0,1,NFFT/2+1);
        freq_res=FS/NFFT;

        Spec = fft(ER,NFFT,2); %size(ER,2); 
        
        mepochs = mean(Spec,1)/(NFFT-1); 
        
        %% Bootstrap
        if opts.bootstrap
            tic;
            L = size(ER,2);
            w = tukeywin(L,0.02/max(time_vec_rec))'; % use a window at the edges to prevent spectral leakage
            Spec = fft(w.*ER,NFFT,2);
            
            strapboot=EFR_straptheboot3(1,Spec,name);
            if opts.save_bootstrapped
                bsname=[savedir cell2mat(condition(condi))];
                save([bsname,'_bootstrap',num2str(channo),'.mat'],'strapboot')
            end
            
            boot_resp_sig=strapboot.bootspecs;
            boot_resp_noise=strapboot.bootnoise;
            boot_resp=(mean(abs(boot_resp_sig),1)-mean(abs(boot_resp_noise),1));
            boot_resp_all = abs(boot_resp_sig) - mean(abs(boot_resp_noise),1);
            boot_resp_all = boot_resp_all./(NFFT-1);
            boot_resp_noise_all = abs(boot_resp_noise)./(NFFT-1);
            
            boot_phase=angle(mean(Spec,1)); %mean(angle(Spec),1);
            mepochs_bs=boot_resp./((NFFT-1));
            epochs_waveform_bs(condi,chani,:) = ifft(boot_resp.*exp(1j*boot_phase),NFFT,'symmetric');

            fprintf("Bootstrapping took %f secs\n",toc);
        end

        %% saving
        epochs_waveform(condi,chani,:) = mean(ER,1);
        epochs_waveform_pos(condi,chani,:) = mean(ERp,1);
        epochs_waveform_neg(condi,chani,:) = mean(ERn,1);
        epochs_all(condi,chani,:)=mepochs;
        if opts.bootstrap
            mepochs_all_bs(condi,chani,:)=mepochs_bs;
            epochs_all_bs(condi,chani,:,:)=boot_resp_all;
            epochs_noise_all_bs(condi,chani,:,:)=boot_resp_noise_all;
        end
        data_all{condi}(chani,:)=Data;
        
        threshold_all{condi,chani} = Thresh;
        nepochs_all(condi,chani,:) = size(ER,1);
    end   
end

czi=find(channels==48);
mepochs_all=2*abs(epochs_all(:,czi,1:NFFT/2+1));
mepochs_all=mean(mepochs_all,2);
mepochs_all=permute(mepochs_all,[1,3,2]);

mepochs_waveform=mean(epochs_waveform(:,czi,:),2);
mepochs_waveform=permute(mepochs_waveform,[1,3,2]);
mepochs_waveform_pos=mean(epochs_waveform_pos(:,czi,:),2);
mepochs_waveform_pos=permute(mepochs_waveform_pos,[1,3,2]);
mepochs_waveform_neg=mean(epochs_waveform_neg(:,czi,:),2);
mepochs_waveform_neg=permute(mepochs_waveform_neg,[1,3,2]);

if opts.bootstrap
    mepochs_all_bs=2*abs(mepochs_all_bs(:,czi,1:NFFT/2+1));
    mepochs_all_bs=mean(mepochs_all_bs,2);
    mepochs_all_bs=permute(mepochs_all_bs,[1,3,2]);

    mepochs_waveform_bs=mean(epochs_waveform_bs(:,czi,:),2);
    mepochs_waveform_bs=permute(mepochs_waveform_bs,[1,3,2]);

    epochs_all_bs=2*abs(epochs_all_bs(:,czi,:,1:NFFT/2+1));
    epochs_all_bs=mean(epochs_all_bs,2);
    epochs_all_bs=permute(epochs_all_bs,[1,3,4,2]); % cond, reps, samples
    
    epochs_noise_all_bs=2*abs(epochs_noise_all_bs(:,czi,:,1:NFFT/2+1));
    epochs_noise_all_bs=mean(epochs_noise_all_bs,2);
    epochs_noise_all_bs=permute(epochs_noise_all_bs,[1,3,4,2]); % cond, reps, samples
end

threshold = {threshold_all{:,czi}};
nepochs = nepochs_all(:,czi,:);
nepochs=permute(nepochs,[1,3,2]);

%% sum amplitudes of harmonics
harmonics=4;

% Use the bootstrapped signal
if opts.bootstrap
    amplitudes_bs=zeros(length(condition),size(epochs_all_bs,2));
    amplitudes_bs_1std=zeros(length(condition),size(epochs_all_bs,2));
    amplitudes_bs_mean=zeros(length(condition),1);
    amplitudes_bs_std=zeros(length(condition),1);
    amplitudes_bs_1std_mean=zeros(length(condition),1);
    amplitudes_bs_1std_std=zeros(length(condition),1);
    for i=1:length(condition)
        if length(f0)==length(condition)
            f0t=f0(i);
        else
            f0t=f0;
        end

        for repi=1:size(epochs_all_bs,2)
            for j=1:harmonics
                peak(j)=fix(j*f0t*NFFT/FS+1);
                amp_dif=epochs_all_bs(i,repi,peak(j));
                if amp_dif>0
                    amplitudes_bs(i,repi)=amplitudes_bs(i,repi)+amp_dif;
                end
                if amp_dif>std(epochs_noise_all_bs,0,3) % keep only the peaks above 1 std
                    amplitudes_bs_1std(i,repi)=amplitudes_bs_1std(i,repi)+amp_dif;
                end
            end
        end
    end
    amplitudes_bs_mean = mean(amplitudes_bs,2); % mean amplitudes
    amplitudes_bs_std = std(amplitudes_bs,0,2); % std of amplitudes
    amplitudes_bs_1std_mean = mean(amplitudes_bs_1std,2); % mean amplitudes (considering peaks above 1 std)
    amplitudes_bs_1std_std = std(amplitudes_bs_1std,0,2); % std of amplitudes (considering peaks above 1 std)
end

% Do the same from the original signal
peak_res=1 * freq_res;     %Hz
floor_res=5 * freq_res;    %Hz
gap_res=1 * freq_res;      %Hz

peak_bins=fix((peak_res/freq_res)/2);
floor_bins=fix((floor_res/freq_res)/2);
gap=fix((gap_res/freq_res)/2);

amplitudes=zeros(length(condition),1);
for i=1:length(condition)
    if length(f0)==length(condition)
        f0t=f0(i);
    else
        f0t=f0;
    end
    
    for j=1:harmonics
        peak(j)=fix(j*f0t*NFFT/FS+1);
        
        peaki(j,:)=peak(j)-peak_bins:peak(j)+peak_bins;
        floori(j,:)=[peak(j)-peak_bins-gap-floor_bins:peak(j)-peak_bins-gap ...
            peak(j)+peak_bins+gap:peak(j)+peak_bins+gap+floor_bins];
        amp_dif=max(mepochs_all(i,peaki(j,:)))-mean(mepochs_all(i,floori(j,:)));
        if amp_dif>0
            amplitudes(i)=amplitudes(i)+amp_dif;
        end
    end
end

%% make output structure
output.freq=freq;
output.time=time_vec_rec;
output.nepochs=nepochs;
output.threshold=threshold;
output.freq_sig=mepochs_all;
output.time_sig=mepochs_waveform;
output.time_sig_pos=mepochs_waveform_pos;
output.time_sig_neg=mepochs_waveform_neg;
if opts.bootstrap
    output.freq_sig_boot=mepochs_all_bs;
    output.time_sig_boot=mepochs_waveform_bs;
    
    output.EFRs_boot=amplitudes_bs_mean;
    output.EFRs_boot_std=amplitudes_bs_std;
    output.EFRs_boot_1std=amplitudes_bs_1std_mean;
    output.EFRs_boot_1std_std=amplitudes_bs_1std_std;
end
output.EFRs=amplitudes;

%% plot
if opts.plot
    set(0,'DefaultTextInterpreter','none');
    for i=1:length(condition)
        if length(f0)==length(condition)
            f0t=f0(i);
        else
            f0t=f0;
        end
        figure
        plot(freq,mepochs_all(i,:),'b','linewidth',1.2,'DisplayName',condition{i})
        for j=f0t:f0t:4*f0t
            hold on
            [~,ji]=min(abs(freq-j));
            plot(j,mepochs_all(i,ji),'r*','linewidth',1.2)
        end
        set(gca,'Fontsize',20)
        xlim([0,610])
        ylim([0, 0.05])
        xlabel('Frequency [Hz]')
        ylabel('EFR amplitude ')
        if strcmp(condition{i},'processed_1000')
            title(['Condition ',num2str(i),' - ','CS_10,0,0',' (Channel=Cz (', num2str(channels(czi)),'))'])
        elseif strcmp(condition{i},'processed_1300')
            title(['Condition ',num2str(i),' - ','CS_13,0,0',' (Channel=Cz (', num2str(channels(czi)),'))'])
        elseif strcmp(condition{i},'processed_700')
            title(['Condition ',num2str(i),' - ','CS_7,0,0',' (Channel=Cz (', num2str(channels(czi)),'))'])
        else
            title(['Condition ',num2str(i),' - ',condition{i},' (Channel=Cz (', num2str(channels(czi)),'))'])
        end
        grid on
        
        h=gcf;
        set(h,'PaperPositionMode','auto');
        set(h,'PaperOrientation','landscape');
        set(h,'Position',[50 50 1200 400]);
        print(gcf, '-dpng', '-r300',[savedir condition{i} '.png'])
    end
end
