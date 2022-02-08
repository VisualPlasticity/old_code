
%function ca_analysis;
% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all;global info;

%% load calcium signal and open recorded image
[f,p]=uigetfile('.signals','load calcium signal data .signals');
cd(p);
load(f,'-mat');
if strfind(f,'pick')
    load([f(1:strfind(f,'pick')-1) '.signals'], '-mat');
end
if strfind(f,'_memmap')
    fname=f(1:strfind(f,'_memmap')-1)
end
if strfind(f,'x')
    temp=uigetfile('*.sbx','Pick the correspondent recording file');
    fname=strtok(temp,'.');
end

if ~exist('sig')
    sig=sig_chunk;
end
if size(sig,1)<size(sig,2)
    sig=sig';
end
total=size(sig,2);

if ~exist('Cor')
    Cor=1:total;
else
    bad=setdiff(1:total,Cor);
end

sig=sig(:,Cor);
%% sync with stimulus data and plot

movement=preprocess(fname);
CAframeHz =info.resfreq/info.recordsPerBuffer;
CAlineHz = info.config.lines;
stim=zeros(1,info.max_idx); %just for plotting purpose,assign stimtype to all the ON frames
for i=1:(numel(info.frame)/2)
    stim(info.frame(i*2-1):info.frame(i*2))=info.stimtype(i);
end

if movement
    try
        load([fname '.ball'],'-mat');
    catch
        load([fname '_ball'],'-mat');
        speed=abs(ball(1:end-1))/192*2./diff(time)'; %size of the imaging area 192pixel, 2cm
    end
    % load([fname '.eye'],'-mat');
    sp=conv(speed,ones(1,ceil(CAframeHz*1))/ceil(CAframeHz*1),'same'); %conv the speed into with a 1s filter
    velocity=interp1(0:1/2:info.max_idx, sp(1:info.max_idx*2+1),1:info.max_idx); % downsampling the speed into the frame num
else
    velocity=[];
end
hsig=sigplt(sig,[stim/max(stim);velocity/max(velocity)],Cor); %
% baseline=prctile(sig,20,1);
% sig=sig./repmat(baseline,size(sig,1),1)-1;

%% set paramenter for stimulus
ncell=size(sig,2);
prestim=ceil(CAframeHz*1); %use 1s baseline
stimON=info.frame(2)-info.frame(1); % StimON duration
seg=prestim+min(info.frame(3:2:end)-info.frame(1:2:(end-2))); % segment=prestim+info.frame(TTL ON+TTL off )

Var=numel(unique(info.stimtype));% calculate types of stimulus, orientations, contrast, etc
rep=floor(min(2*numel(info.stimtype),numel(info.frame))/2/Var); %calculate repetitions

sigT=zeros(seg*rep,Var,ncell); % build sigT based on signals from sig and organized according to segment
frame_on=zeros(Var,rep); %tag frame_on of each seg
window=prestim:prestim+stimON;
%% reorganize stimulus according to stimulus type
last=info.frame(end-1)+seg-prestim;
if last>size(sig,1)
    sig(end+1:last,:)=0;
    disp(sprintf('padding zeros to the last %d frames', last-size(sig,1)));
end

for j=1:Var
    ori=find(info.stimtype==j);
    on=ori*2-1;
    frame_on(j,:)=info.frame(on(1:rep)); % the frame marking the prestim timepoint size 1*rep
    temp=ones(seg,1)*frame_on(j,:)+ (0:seg-1)'*ones(1,rep)-prestim; % array seg*rep6
    sigT(:,j,:)=sig(temp,:);   %(seg*rep)*Var*ncell
end

sigT=reshape(sigT,seg,rep,Var,ncell);
%% select cells based on their peak responses

if ~isempty(strfind(f,'memmap')) & isempty(strfind(f,'pick2'))
    hcell=openfig([f(1:findstr(f,'cell')+3) '.fig']);
    if ~exist('bad','var')
        bad=pick(f);
    end
    [bad2,hpick]=pick2(sigT,window,Cor);
    bad=unique([bad bad2]);
    if ~isempty(bad)
        Cor = setdiff(1:ncell,bad);
        ncell= numel(Cor);
        newf=[f(1:end-8) 'pick2ed' num2str(ncell)];
        
        axesObjs = get(hcell, 'Children');  %axes handles
        dataObjs = get(axesObjs, 'Children');
        temp=dataObjs(total+1-bad);
        delete(temp);
        savefig(hcell,[newf '.fig']);
        
        sig=sig(:,Cor);
        sigT=sigT(:,:,:,Cor);
        save([newf '.signals'],'Cor');
    end
end

%% syn with ball/eye motion
if movement
    run=[];%frame_on=zeros(Var,rep)
    for k=1:stimON
        run=cat(3,run,velocity(frame_on'+k-1));% run:rep,Var,stimON
    end
    [matrix,hrun]=runplt(run);
    %[matrix,hrun]=runplt2(run);
else
    matrix=[];
end
%% parse sigF according to different types of info.var

if info.steps(2)>1
    [ sigF,variant,hparse]=parseF(sigT,window,matrix);

else
    sigF=sigT; % sigF:seg,rep,Var,ncell
    variant=strcmp(info.var(1),'Orientation');
end
save([fname '_' num2str(numel(Cor)) 'signals.mat'],'sigF','matrix','window','Cor');


%% plot figures and calculate Orientation selection or Contrast invariant
save([fname '_' num2str(numel(Cor)) 'signals.mat'],'sigF','matrix','window','Cor');

drawingoption=(~isempty(matrix))*10+variant;

switch drawingoption
    case 0%'no running with contrast map'
        hsigF=sigFplt(sigF,[],window,Cor);  % sigF:seg,rep,Var,ncell
        [SI.peak,~,SI.error]=cal_ER(sigF,window);%peak: 1*Var*ncell
        
        hpk(1)=fitplt(SI.peak,SI.error,Cor);
        hpk(1).Name='Contrast tuning under best orientation';
        
    case 1 %'no running with orientation map'
        [SI.peak,~,SI.error]=cal_ER(sigF,window);  %peak:1*
        hsigF=sigFplt(sigF,[],window,Cor);  % sigF:seg,rep,Var,ncell
        
        hpk(1)=polarplt(SI.peak,Cor);
        hpk(1).Name='ori tuning under best contrast';
        
        SI.OSI=calOSI(SI.peak);
        SI.DSI=calDSI(SI.peak);
        SI.gOSI=calgOSI(SI.peak);
        
        hpk(2)=figure('Name','orientation tuning under best contrast');hold on;        
        fitplt(SI.peak,[],Cor);
        
    case 2 %'no running with both contrast&orientation map'
        hsigF=sigVarplt(sigF,window,[],Cor);  % sigF:seg,rep,Var,ncell
        [SI.peak,~,SI.error]=cal_ER(sigF,window);%peak:1*Var*ncell
        SI.peak=reshape(SI.peak,info.steps(2),info.steps(1),ncell);
        SI.error=reshape(SI.error,info.steps(2),info.steps(1),ncell);
        
        hpk(1)=fitplt(SI.peak,SI.error,Cor);           %peak: 1*Var*ncell
        hpk(1).Name='Orientation tuning under diff contrast-running';
        
        hpk(2)=polarplt(SI.peak,Cor);
        
        SI.OSI=calOSI(SI.peak);
        SI.DSI=calDSI(SI.peak);
        SI.gOSI=calgOSI(SI.peak);
        
        hpk(3)=fitplt(SI.gOSI,[],Cor);
        hpk(3).Name='gOSI under diff contrast';
        
    case 10  %'running with contrast'
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,window,matrix);  % sigR:seg,1,Var,ncell
        hsigF=sigFplt(sigF,matrix,window,Cor);  % sigF:seg,rep,Var,ncell
        SI.peakR=fixpeak(SI.peakR);
        SI.peakS=fixpeak(SI.peakS);

        
        hpk(1)=figure('Name','Contrast tuning under best orientation-running');hold on;
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
        
        
    case 11  %'running with orientation'
        [hsigF,Gd]=sigFplt(sigF,matrix,window,Cor);  % sigF:seg,rep,Var,ncell
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,window,matrix);  % sigR:seg,1,Var,ncell
        hpk(1)=polarplt(cat(1,SI.peakR,SI.peakS),Cor);
        
        
        SI.OSI_R=calOSI(SI.peakR);        
        SI.OSI_S=calOSI(SI.peakS);        
        SI.gOSI_R=calgOSI(SI.peakR);
        SI.gOSI_S=calgOSI(SI.peakS);
        SI.DSI_R=calDSI(SI.peakR);
        SI.DSI_S=calDSI(SI.peakS);  
        
        hpk(2)=figure('Name','OSIcomparison');hold on;
        scatter(SI.gOSI_S(:),SI.gOSI_R(:),'jitter','on', 'jitterAmount', 0.15);
        scatter(SI.OSI_S(:),SI.OSI_R(:),'g','jitter','on', 'jitterAmount', 0.15);
        plot([0 max(SI.gOSI_R(:))],[0 max(SI.gOSI_R(:))],'--');
        xlabel('gOSI still');ylabel('gOSI running');
        
        hpk(3)=figure('Name','DSIcomparison');
        scatter(SI.DSI_S(:),SI.DSI_R(:),'jitter','on', 'jitterAmount', 0.15);
        hold on;
        plot([0 max(SI.DSI_R(:))],[0 max(SI.DSI_R(:))],'--');
        xlabel('DSI still');ylabel('DSI running');
        
        
        
    case 12%'running with both contrast and orientation'
        hsigF=sigVarplt(sigF,window,matrix,Cor);  % sigF:seg,rep,Var,ncell
        [~,SI.peakR,SI.errorR,~,SI.peakS,SI.errorS]= sigFcmp(sigF,window,matrix);  % calculate peak from prestim:prestim+stimON-1
        SI.peakR=reshape(SI.peakR,info.steps(2),info.steps(1),ncell);
        SI.peakR=fixpeak(SI.peakR);
        SI.errorR=reshape(errorR,info.steps(2),info.steps(1),ncell);
        
        
        SI.peakS=reshape(SI.peakS,info.steps(2),info.steps(1),ncell);
        SI.peakS=fixpeak(SI.peakS);
        SI.errorS=reshape(SI.errorS,info.steps(2),info.steps(1),ncell);
        
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
        hpk(1).Name='Orientation tuning under diff contrast-running';
        hpk(2)=polarplt(cat(1,SI.peakR,SI.peakS),Cor);
        % SI.peakR contrast*ori*ncell
        
        SI.OSI_R=calOSI(SI.peakR);        
        SI.OSI_S=calOSI(SI.peakS);
        
        SI.gOSI_R=calgOSI(SI.peakR);
        SI.gOSI_S=calgOSI(SI.peakS);
        
        hpk(3)=fitplt(cat(1,SI.gOSI_R,SI.gOSI_S),[],Cor);
        hpk(3).Name='gOSI under diff contrast and running state';
        
end

%% save figures and parameters
% prompt = 'Do you want more? Y/N [Y]: ';
% str = input(prompt,'s');
% if isempty(str)
%     str = 'Y';
% end

mark=findstr(fname,'_');
date=fname(mark(1)+1:mark(1)+3);
serial=fname(mark(2)+1:mark(2)+3);
folder=fullfile([ date '-' serial ],['picked' num2str(numel(Cor))],['drawingoption' num2str(drawingoption)]);
if ~exist(folder,'dir')
    mkdir(folder);
end

save(fullfile(folder,'peakSI.mat'),'SI','Gd');

for i=1:numel(hsigF)
    saveas(hsigF(i),fullfile(folder,hsigF(i).Name),'png');
end
for i=1:numel(hpk)
    savefig(hpk(i),fullfile(folder,hpk(i).Name));
    saveas(hpk(i),fullfile(folder,hpk(i).Name),'png');
end

if exist('hrun')
for i=1:numel(hrun)
    saveas(hrun(i),fullfile(folder,hrun(i).Name),'png');
end
end

if exist('hparse')c
    for i=1:numel(hparse)
    saveas(hparse(i),fullfile(folder,hparse(i).Name),'png');
    end
end
if exist('hpick')
    saveas(hpick,fullfile(folder,hpick.Name),'png');
end
