function [buf] = play_buf(y, playDeviceID, playChanList, Fs)
% PLAY_BUF
%    playrec_buf(y, playDeviceID , playChanList,  Fs)
%    sends the stimuli and the triggers to the soundcard.
pageSize = 4096; %256;
pageBufCount = 5;
runMaxSpeed = true;
startPoint = 1;
endPoint = size(y,1);
fileChanCount = size(y,2);
if length(playChanList) ~= fileChanCount 
    error('y must have the same number of columns than number of channel in chanList');
end
if(playrec('isInitialised'))
    if(playrec('getSampleRate')~=Fs)
        fprintf('Changing playrec sample rate from %d to %d\n', playrec('getSampleRate'), Fs);
        playrec('reset');
    elseif(playrec('getPlayDevice')~=playDeviceID)
        fprintf('Changing playrec play device from %d to %d\n', playrec('getPlayDevice'), playDeviceID);
        playrec('reset');      
    elseif(playrec('getPlayMaxChannel')<max(playChanList))
        fprintf('Resetting playrec to configure device to use more output channels\n');
        playrec('reset');
    end
end
if(~playrec('isInitialised'))
    fprintf('Initialising playrec to use sample rate: %d, playDeviceID: %d and recDeviceID: %d\n', Fs, playDeviceID);
    playrec('init', Fs, playDeviceID, -1, max(playChanList));
    pause(0.1);
end   
if(~playrec('isInitialised'))
    error ('Unable to initialise playrec correctly');
end
if(playrec('pause'))
    fprintf('Playrec was paused - clearing all previous pages and unpausing.\n');
    playrec('delPage');
    playrec('pause', 0);
end
playrec('delPage');
buf = [];
pageNumList = repmat(-1, [1 pageBufCount]);
firstTimeThrough = true;
for startSample = startPoint:pageSize:endPoint
    endSample = min(startSample + pageSize - 1, endPoint);    
    yp = y(startSample:endSample,:);
    pageNumList = [pageNumList playrec('play', yp, playChanList)]; %#ok<AGROW>
    
    if(firstTimeThrough)
        playrec('resetSkippedSampleCount');
        firstTimeThrough = false;
    else
        if(playrec('getSkippedSampleCount'))
            fprintf('%d samples skipped!!\n', playrec('getSkippedSampleCount'));
            firstTimeThrough = true;
        end
    end
    if(runMaxSpeed)
        while(playrec('isFinished', pageNumList(1)) == 0) %#ok<UNRCH>
        end
    else
        playrec('block', pageNumList(1));
    end
    lastRecording = playrec('getRec', pageNumList(1));
    if(~isempty(lastRecording))
        buf = [buf; lastRecording]; %#ok<AGROW>
    end
    playrec('delPage', pageNumList(1));
    pageNumList = pageNumList(2:end);
end
while ~isempty(pageNumList>1)
    if(runMaxSpeed)
        while(playrec('isFinished', pageNumList(1)) == 0) %#ok<UNRCH>
        end
    else
        playrec('block', pageNumList(1));
    end
    
    playrec('delPage', pageNumList(1));
    pageNumList = pageNumList(2:end);
end
playrec('delPage');
end