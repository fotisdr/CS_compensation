function [xp,En,g]=g_70dB(x,cs,win_len,height,prominence,distance)
% The g_70dB function processes an input stimulus based on the desired 
% cochlear-synaptopathy (CS) profile parameters. The envelope processing 
% of a 70 dB-SPL stimulus is applied here, irrespectively of the amplitude 
% of the signal.
% 
% x is the input stimulus - can also be two-dimensional for processing 
% multiple signals, with the second dimension  corresponding to the time 
% dimension
% cs is the CS profile for which the processing will be applied - 3 CS 
% profiles are included: 13,0,0, 10,0,0 and 7,0,0 
% 
% The remaining parameters are optional and are used for the envelope 
% estimation of the input stimulus. If no parameters are given then the 
% Hilbert envelope of the signal is used, otherwise an RMS-based estimation
% of the envelope is performed. 
% win_len is the length of the running RMS window 
% height is the minimum height of envelope peaks to be considered
% prominence is the minimum difference from adjacent peaks
% distance is the minimum distance between peaks

% define the CS parameters of the non-linear processing function g
a = 1;
Emax_70dB = 0.145774430020356; % maximum amplitude of the envelope at 70 dB SPL
if length(cs) == 1
    if cs == 1300
        b = 50.45;
        c = 4.8476;
    elseif cs == 1000
        b = 270.1903;
        c = 30.6769;
    elseif cs == 700
        b = 1999.1160;
        c = 260.1633;
    else
        error('CS profile not supported')
    end
elseif length(cs) == 3
    if cs(1) == 13 & cs(2) == 0 & cs(3) == 0
        b = 50.45;
        c = 4.8476;
    elseif cs(1) == 10 & cs(2) == 0 & cs(3) == 0
        b = 270.1903;
        c = 30.6769;
    elseif cs(1) == 7 & cs(2) == 0 & cs(3) == 0
        b = 1999.1160;
        c = 260.1633;
    else
        error('CS profile not supported')
    end
else
    error('CS profile not supported')
end

if nargin<3
    win_len=0;
end

% the processing is applied along the second dimension
if iscolumn(x)
    x=x';
    inv_flag=1;
else
    inv_flag=0;
end

E=zeros(size(x));
if win_len
    for n=1:size(x,1)
        % compute the envelope of the signal using the RMS-based method
        E(n,:)=abs(envelope(x(n,:),win_len,'rms'));
        if exist('smoothdata')
            E(n,:)=smoothdata(E(n,:),'movmean',win_len);
        else
            E(n,:)=smooth(E(n,:),'moving',win_len);
        end
    end
    
    % find the local maxima (peaks) of the envelope
    peaks=1;
    minima=[];
    cntr=1;
    for i=3:size(E,2)
        if E(n,i-1)>height & E(n,i-1)>E(n,i) & ...
                E(n,i-1)>E(n,i-2) & i-1>=peaks(cntr)+distance
            [~,peak]=min(E(n,peaks(cntr):i-1));
            if E(n,i-1)>prominence*E(n,peaks(cntr)+peak-1)
                minima=[minima peaks(cntr)+peak-1];
                cntr=cntr+1;
                peaks(cntr)=i-1;
            end
        end
    end
else
    % compute the Hilbert envelope of the signal if no parameters are given
    for n=1:size(x,1)
        E(n,:)=abs(envelope(x(n,:)));
    end
end

% compute the local minima (troughs) of the envelope
minima(1)=1;
if minima(end)~=size(E,2)
    minima=[minima,size(E,2)];
end
% if win_len is not given then the whole envelope is processed at once
if win_len
    peaks=peaks(2:end);
    minimum=1;
    for i=1:length(peaks)
        if E(n,peaks(i))>=prominence*E(n,minima(i+1))
            minimum=[minimum minima(i+1)];
        end
    end
    if minimum(end)~=size(E,2)
        minimum=[minimum size(E,2)];
    end
    minima=minimum;
end

%% CS compensation
if any(cs~=[13,3,3])
    % compute the processing function between each pair of envelope minima
    for j=1:length(minima)-1
        index=minima(j):minima(j+1);
        x_segment=x(:,index);
        
        for n=1:size(x_segment,1)
            Emax(n)=max(abs(x_segment(n,:))); % temporal max of the envelope
            Emaxr(n)=Emax(n)./Emax_70dB; % relative max to the 70-dB-SPL max used in the optimization
            if Emaxr(n) == 0
                Emaxr(n) = 1e-20;
            end

            % normalize the envelope to the local maximum
            En(n,index)=E(n,index).*Emax(n)./max(abs(E(n,index)));
            % recalibrate the "a" parameter to ensure that the peak stays exactly the same
            at(n)=a./max(1./(1+exp(-b.*En(n,index)./Emaxr(n)+c)));
            % compute the processing function g over the stimulus segment
            g(n,index)=at(n)./(1+exp(-b.*En(n,index)./Emaxr(n)+c));
        end
    end
    % remove nan or inf values (due to division by zero)
    for n=1:size(x_segment,1)
        g(n,:)=fillmissing(g(n,:),'linear');
    end
    % process the stimulus using the computed processing function g
    xp=g.*x;
else
    xp=x; % if 13,3,3 is given apply no processing
end

if inv_flag
    xp=xp';
end
