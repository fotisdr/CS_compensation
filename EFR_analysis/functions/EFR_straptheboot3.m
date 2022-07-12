function [ strapboot ] = EFR_straptheboot3(art,EFR,fname) % fname|string, nh_efr|cell array including strucutres, subj|integer, subjn|string, cond|string, badchan|vector
% Bootstrapping method from the Zhu et al. ASA publication
% Sarineh Keshishzadeh, UGent
%%
%%*************************************************************************

pars = parsconfig_efr('bootstrap');
snum = sum(double(fname).^2); % for the random number generator; turns data name into integer
data(1,:,:)=EFR;
NFFT = 2^nextpow2(size(data,3));
% fprintf('\n:: %d epochs were extracted\n', size(data,2))

%%  Start loop for all channels in pars.nchannels

for n = pars.nchannels  %for each channel!! calculate the noise floor for each channel!
%     tic;
%     fprintf('\n:: Compute bootspecs for channel %d\n',n)
    
%     bootspecs = zeros(pars.nDraws,NFFT/2+1); % final size of matrix, all zeros
    Freq = pars.fs/2*linspace(0,1,NFFT/2+1); % # Frequencies for FFT results
    
    for j = 1:pars.nDraws
        
        %fprintf('\n:: Draw %d out of %d',j,pars.nDraws)
        rng(1+n+j+snum); % seed for the random number generator: 1 + channel + draw + id number
        draws = randi(size(data,2),1,size(data,2));
        draws = squeeze(data(n,draws,:));
%         efr = zeros(size(data,2),NFFT);
        
      
            efr=mean(draws,1);  % compute the fft for each epoch
   
%         MDFT=2*abs(fft(efr,NFFT)/size(draws,2)); % average the ffts over all epochs
%         bootspecs(j,:) = MDFT(1:NFFT/2+1);
        bootspecs(j,:)=efr;
    end
%     toc
    
    %% Noisefloor (Phase flip method, faster)
    if art == 1
%     tic;
%     fprintf('\n:: Using Phase flip method')    
%     fprintf('\n:: Compute bootnoise for channel %d\n',n)
    
%     bootnoise = zeros(pars.nDraws_noise,NFFT/2+1);
    
    for j = 1:pars.nDraws_noise
        
        %fprintf('\n:: Draw %d out of %d',j,pars.nDraws_noise)
        
        rng(2+n+j+snum); % seed for the random number generator: 2 + channel + draw * subj id
        draws_one = randi(size(data,2),1,floor(size(data,2)/2));
        rng(3+n+j+snum); % seed for the random number generator: 3 + channel + draw * subj id
        draws_two = randi(size(data,2),1,floor(size(data,2)/2));
        
        draws_one = squeeze(data(n,draws_one,:));
        draws_two = squeeze(data(n,draws_two,:))*-1;
        draws = cat(1,draws_one,draws_two);
        
        
        
       efr1=mean(draws,1);   % compute the fft for each epoch
        
%         MDFT=2*abs(fft(efr,NFFT)/size(draws,2));% average the ffts over all epochs
%         bootnoise(j,:) = MDFT(1:NFFT/2+1);
 bootnoise(j,:) =efr1;
        
    end
    end
%     fprintf('\n')
%     toc
    %% Noisefloor - Phase randomization (takes much longer, similar results)
    
    if art == 2
%     tic;
%     fprintf('\n:: Using Phase randomization method')
%     fprintf('\n:: Compute bootnoise for channel %d\n',n)
    bootnoise = zeros(pars.nDraws_noise,NFFT/2+1);
    
    for j = 1:pars.nDraws_noise
        
        %fprintf('\n:: Draw %d out of %d',j,pars.nDraws_noise)
        
        rng(4+n+j+snum); % seed for the random number generator: 4 + channel + draw * subjnumber
        draws = randi(size(data,2),1,size(data,2));
        draws = squeeze(data(n,draws,:));
        DFT = zeros(size(data,2),NFFT);
        
        for k = 1:size(draws,1)
            
            
            DFT(k,:)=fft(draws(k,:),NFFT)/size(draws,2);   % compute the fft for each epoch
            R = abs(DFT(k,:));
            theta = 2*pi.*rand(1,size(R,2));
            %and the statement
            DFT(k,:) = R.*exp(i*theta);
            
        end
        MDFT=mean(DFT,1); % average the ffts over all epochs
        bootnoise(j,:) = 2*abs(MDFT(1:NFFT/2+1));
        
    end
    end
%     fprintf('\n')
%     toc
    %%
    if length(pars.nchannels) == 1
        strapboot.bootspecs = bootspecs;
        strapboot.bootnoise = bootnoise;
    else
        strapboot.bootspecs{1,n} = bootspecs;
        strapboot.bootnoise{1,n} = bootnoise;
    end
    strapboot.freq = Freq;
    
end

%beep
% fprintf('\n:: Function has finished\n')
end
