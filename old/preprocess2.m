function movement=preprocess2(fname)
global info;

load([fname '.align'],'m','-mat');
sbxread(fname,1,1);
if info.frame(1)==0
    info.frame=info.frame(2:end);
    info.line=info.line(2:end);
    %    save([fname '.mat'],'info','-append');
end
%% match stim
if ~isfield(info,'stimtype')   %read psychostim stim info to match TTL event recorded in
    file=dir([fname '.sbx']);
    recordedtime=file.date;
    [day,time]=strtok(recordedtime)
    %[f2,p2]=uigetfile(['\\mps-pc53\2pdata\stimulus\' day],['find stim info just before' day time]); %'
    [f2,p2]=uigetfile(['\\mps-pc38\Documents\Matlab\stimulus\' day],['find stim info just before' time]);%'
    load(fullfile(p2,f2), '-mat');
    info.sttype = intanSyncData(:,1);% the variable IntanSyncData has all of the useful information
    info.counting= intanSyncData(:,2);
    info.time = intanSyncData(:,3);
	info.stimtype=info.sttype(info.counting==1);%find(diff(sttype),1) can define the length of ON and OFF, so, info.stimtype reads the first frame of each stimulus
%     info.stimtype=sttype(1: find(diff(sttype),1):end);
    %find(diff(sttype),1) can define the length of ON and OFF, info.stimtype reads the first frame of each stimulus
    %info.stimtype=sttype(logical([1;diff(sttype)]));
    info.stimtype = max(info.stimtype) -info.stimtype;

    
    
    info.var=Var1Str{Var1Val};
    info.steps=nSteps1;
    if Var2Val>1
        info.var={info.var,Var2Str{Var2Val}};
        info.steps=[info.steps,nSteps2];
    else
        info.var={info.var,'none'};
        info.steps=[info.steps,1];
        %        save([fname '.mat'],'info.stimtype','info.sttype','info.var','-append');
        
    end
          
    
    short=abs(2*numel(info.stimtype)-numel(info.frame));
    if short
        figure;plot(info.time,info.counting/max(info.counting),'o--',...
        info.time,info.sttype/max(info.sttype),'*-');
    disp(info.stimtype(1:5))
    disp(info.frame(1:10))
        info.stimtype = max(info.stimtype) -info.stimtype;

        xlim([0 30]);
        warning=sprintf('stimtype# %d and recorded# %d do not match,check if necessary',2*numel(info.stimtype),numel(info.frame));
        disp(warning);
    else
       % save([fname '.mat'],'info','-append');
        disp('stim and recording matched and saved!');
    end
end

if exist([fname '_ball.mat'])
    movement=1;
    if ~exist([fname '.ball'])
        sbxballmotion([fname '_ball.mat']) ;
    end
    display('with ball info, will analyze running later');
else
    display('no ball info');
    movement=0;
end
%% Do not analysis eye info at the moment
% if exist([fname '_eye.mat'])
%     movement=movement+10;
%     if ~exist([fname '.eye'])
%         sbxeyemotion([fname '_eye.mat']) ;
%     end
%     display('with movement info, will analyze running later');
% else
%     display('no eye info');
% end