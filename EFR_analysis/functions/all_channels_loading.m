function all_channels_loading(condition,Chanselect,CR1,CR2,datadir,savedir,chunks)
% Load the data of the selected channels (defined in Chanselect) for a file
% (defined in condition). The data are normalised by the two reference 
% channels CR1 and CR2 and then saved under savedir.

if nargin < 5
    datadir='';
else
    datadir = [datadir filesep];
end
if nargin < 6
    savedir='';
else
    savedir = [savedir filesep];
end
if nargin < 7
    chunks=100;
end

for name_condition=1:length(condition)
    name=[datadir cell2mat(condition(name_condition))];
    savename=[savedir cell2mat(condition(name_condition))];
    
    eval(['load ',num2str(name),'_',num2str(1),'_',num2str(chunks),'.mat']);
    ce = data.Head.NRec;

    %65 & 66 are the EX1 and Ex2 channels respectively
    %if you want to average across Ex1 and Ex2: CR1=65 and CR2=66
    %if you want to use Ex1 alone: CR1=65 and CR2=65
    %if you want to use Ex1 alone: CR1=65 and CR2=65
    
    %for Channo=1:64 %if you want to run all channels=> will take a long time!!
    for Channo=Chanselect
        
        %make a vector that has the same length as the total file
        %store channel concatenated and referenced on mean 1 and 2 of exernal channels;
        %channel 33 and 34 are the reference electrodes
        display(['Processing channel ',num2str(Channo),'/70'])
        %concatenate and reference the data of one channel
        for n=0:ceil(ce/chunks)-1
            if n==ceil(ce/chunks)-1
                eval(['load ',num2str(name),'_',num2str(1+n*chunks),'_',num2str(ce),'.mat']);
                Lblsize=size(data.Record,2);
                Ch(1,1+n*blsize:n*blsize+size(data.Record,2))=data.Record(Channo,:)-(data.Record(CR1,:)+data.Record(CR2,:))/2;
                display(['concatenate no ',num2str(n),'/',num2str(ceil(ce/chunks)-1)])
                clear data
            elseif n==0
                eval(['load ',num2str(name),'_',num2str(1+n*chunks),'_',num2str((n+1)*chunks),'.mat']);
                Ch=zeros(1,size(data.Record,2)*ceil(ce/chunks));
                Ch(1,1+n*size(data.Record,2):(n+1)*size(data.Record,2))=data.Record(Channo,:)-(data.Record(CR1,:)+data.Record(CR2,:))/2;
                blsize=size(data.Record,2);
                display(['concatenate no ',num2str(n),'/',num2str(ceil(ce/chunks)-1)])
                clear data
            else
                eval(['load ',num2str(name),'_',num2str(1+n*chunks),'_',num2str((n+1)*chunks),'.mat']);
                Ch(1,1+n*size(data.Record,2):(n+1)*size(data.Record,2))=data.Record(Channo,:)-(data.Record(CR1,:)+data.Record(CR2,:))/2;
                display(['concatenate no ',num2str(n),'/',num2str(ceil(ce/chunks)-1)])
                clear data
            end
        end
        Ch(n*blsize+Lblsize+1:end)=[];
        
        save([savename,'Ch',num2str(Channo),'.mat'],'Ch');
        clear Ch
    end
end

