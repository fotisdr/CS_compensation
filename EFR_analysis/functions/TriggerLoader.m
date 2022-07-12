function TriggerLoader(name,savename,condition,chunks)
% Load triggers from the mat files extracted from a recorded bdf file 
% (name) and save them under the savename location.

if nargin < 2
    savename=name; 
end
if nargin < 3
    condition=name; 
end
if nargin < 4
    chunks=100;
end
eval(['load ',num2str(name),'_',num2str(1),'_',num2str(chunks),'.mat']);
ce = data.Head.NRec;

Channo=1; % channel to load the triggers from

%make a vector that has the same length as the total file
%concatenate the data of the trigger channel
for n=0:ceil(ce/chunks)-1
    if n==ceil(ce/chunks)-1 
        eval(['load ',num2str(name),'_',num2str(1+n*chunks),'_',num2str(ce),'.mat']); 
        Lblsize=size(data.Record,2);
        Trigs(1,1+n*blsize:n*blsize+size(data.Triggers,2))=data.Triggers(Channo,:);
        display(['concatenate no ',num2str(n),'/',num2str(ceil(ce/chunks)-1)])
        clear data
    elseif n==0
        eval(['load ',num2str(name),'_',num2str(1+n*chunks),'_',num2str((n+1)*chunks),'.mat']);
        Trigs=zeros(1,size(data.Record,2)*ceil(ce/chunks));
        Trigs(1,1+n*size(data.Record,2):(n+1)*size(data.Record,2))=data.Triggers(Channo,:);
        blsize=size(data.Record,2);
        display(['concatenate no ',num2str(n),'/',num2str(ceil(ce/chunks)-1)])
        clear data
    else
        eval(['load ',num2str(name),'_',num2str(1+n*chunks),'_',num2str((n+1)*chunks),'.mat']);
        Trigs(1,1+n*size(data.Record,2):(n+1)*size(data.Record,2))=data.Triggers(Channo,:);
        display(['concatenate no ',num2str(n),'/',num2str(ceil(ce/chunks)-1)])
        clear data
    end
end

Trigs(n*blsize+Lblsize+1:end)=[];
save([savename 'Trigs',condition,'.mat'],'Trigs');

end