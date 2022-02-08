% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
%% load calcium signal and plot
[f,p]=uigetfile('.signals','load calcium signal data .signals');
cd (p);load(f);
load(f,'-mat')

figure;imshow(m)
 if(exist([a '.segment'],'file'))
        load([a '.segment'],'vert','-mat');
        ch = get(gca,'children');
        for(i=1:length(vert))
            if ~isempty(vert{i})
                temp=vert{i}';
                text(temp(1),temp(2),num2str(i),'Fontsize',12);
                n = n + 1;
                drawnow
            end
        end
    end
%% load stimulus data and sync
fname=f(1:(findstr(f,'_memmap')-1));
sbxread(fname,1,1);
global info;
CAframeHz =info.resfreq/info.recordsPerBuffer;
CAlineHz = info.config.lines;


if info.frame(1)==0
    info.frame=info.frame(2:end);
    info.line=info.line(2:end);
    save([fname '.mat'],'info','-append');
end

if ~isfield(info,'stimtype')   %read psychostim stim info to match TTL event recorded in
    file=dir([fname '.sbx']);
    recordedtime=file.date;
    [day,time]=strtok(recordedtime)
    [f2,p2]=uigetfile(['..\..\stim\' day],['find stim info just before' time]);
    load(fullfile(p2,f2), '-mat');
    sttype = intanSyncData(:,1);% the variable IntanSyncData has all of the useful information
    info.stimtype=sttype(1: find(diff(sttype),1):end);%find(diff(sttype),1) can define the length of ON and OFF, so, info.stimtype reads the first frame of each stimulus
end

if (2*numel(info.stimtype))==numel(info.frame)
    save([fname '.mat'],'info','-append');short=0;
else
    warning=sprintf(['stimtype#%d and recorded#%d do not match,check if necessary',numel(info.stimtype),numel(info.frame)]);
    disp(warning);
    short=
end
%info.stimtype=ceil(info.stimtype/4);     % mod(info.stimtype-1,8)+1;


%% set paramenter for stimulus
ncell=size(C_df,1);
prestim=ceil(CAframeHz*1); %use 1s baseline

orientations=numel(unique(info.stimtype));% calculate types of orientations
rep=numel(info.stimtype)/orientations-short; %calculate repetitions 
stimON=info.frame(2)-info.frame(1); % StimON duration
sigT=zeros(orientations,seg,ncell); % build sigT based on signals from sig and organized according to segment
frame_on=zeros(orientations,rep); %tag frame_on of each seg
peak=zeros(orientations,ncell); %  

%% 

stim=zeros(1,info.max_idx); %just for plotting purpose,assign stimtype to all the ON frames 
for i=1:(numel(info.frame)/2)
    stim(info.frame(i*2-1):info.frame(i*2))=info.stimtype(i);
end

fn=fname;
fn(strfind(fname,'_'))='-';

figure('Name',[fn 'data with stim &ball movement'], 'position',[200 100 2400 1000]);
n=ncell+movement;
margin=0.05;
xpos=margin;
ypos=((1:n)-1)/n;
xwidth=1-margin*2;
yheight= 1/n-margin*2/n;
ymax=max(sigF(:));
ymin=min(sigF(:));

subplot('position',[xpos(1) ypos(1) xwidth yheight]);
plot([ 1;prestim],[orientations;orientations],'-k',[1;1],[0;orientations],'-k','Linewidth',2)
text(0,orientations/2, num2str(orientations), 'HorizontalAlignment','left')
text((1+prestim)/2,orientations-1, num2str(prestim-1), 'HorizontalAlignment','center')

for k=1:ncell
    subplot('position',[xpos ypos(k) xwidth yheight]);hold on; 
    plot(sigF(:,k),'linewidth',2)
    plot(stim,'k','linewidth',1);
    axis([1 info.max_idx ymin max(ymax,orientations)]);
    text(1,orientations/2,['cell' num2str(k)],'HorizontalAlignment','right');
    %axis off
end


if movement
    subplot('position',[xpos ypos(k+1) xwidth yheight]);
    hold on;
 %   plot(velocity);
    plot(stim,'k','linewidth',1);
    axis tight;
    plot( [0;0],[0;10],'-k','Linewidth',2);
    text(-1,5, '10cm/s')
    axis off
end

for j=1:orientations
    ori=find(info.stimtype==j);
    on=ori*2-1;
    frame_on(j,:)=info.frame(on(1:rep));
%     for k=1:rep
%         run(j,k)=mean(velocity(frame_on(j,k):frame_on(j,k)+stimON-1));
%     end
end

%%
figure('Name',[fn 'reorganized'],'position',([ 200 100 1600 1200]));
margin=0.05;
xrow=orientations;
ycol=ncell;
xpos=margin+ ((1:xrow)-1)/xrow* (1-margin*2);
xwidth=(1-margin*2)/xrow* (1-margin*2);
ypos=margin+ ((1:ycol)-1)/ycol* (1-margin*2);
yheight= (1-margin*2)/ycol* (1-margin*2);
ymax=2;
ymin=0;


subplot('position',[xpos(1) ypos(1) xwidth yheight]);
plot([ 0;100],[0;0],'-k',[1;1],[0;2],'-k','Linewidth',2)
text(0,2/2, '2', 'HorizontalAlignment','right')
text(100/2,-1,'100', 'HorizontalAlignment','center')

for i=1:ncell
%    m=vert{i};
%     try
%         line=max(m(:,2)); %line of the scanned cell
        for j=1:orientations
        %    frame_on(j,:)=frame_on(j,:)+(info.line(on)>=line);    % If info.line(on) is larger, means the cell is already scanned before stimON, use the next frame;
            subplot('position',[xpos(j) ypos(i) xwidth yheight]);
            hold on;
            temp=[];
            run_rep=0;
            for k=1:rep
                 baseline=min(sigF((frame_on(j,k)-prestim):frame_on(j,k),i));%plot 2s prestim and use the mim for baseline
                 temp=sigF((frame_on(j,k)-prestim):(frame_on(j,k)-prestim+seg-1),i)-baseline; %prestim +seg 
                plot(temp,'k-');
                sigT(j,:,i)=sigT(j,:,i)+temp';
 %               run(j,k)=mean(velocity(frame_on(k):frame_on(k)+stimON-1));
%                 if run(j,k)> .5
%                     sigT_run(j,:,i)=sigT_run(j,:,i)+temp';
%                     run_rep=run_rep+1;
%                     plot(temp,'r-');
%                 else
%                     sigT_still(j,:,i)=sigT_still(j,:,i)+temp';
%                     plot(temp,'b-');
%                 end
            end
            sigT(j,:,i)=sigT(j,:,i)/k; % count into the new segment and do the average
%             if run_rep>0
%                 sigT_run(j,:,i)=sigT_run(j,:,i)/run_rep;
%             end
%             sigT_still(j,:,i)=sigT_still(j,:,i)/(rep-run_rep);
            
            peak(j,i)=max(sigT(j,prestim+1:end,i));   %peak is calculated from the 15ms from the stimulus onset towards the end, but each frame is 60ms...
%             peak_run(j,i)=max(sigT_run(j,prestim+1:end,i)); 
%             peak_still(j,i)=max(sigT_still(j,prestim+1:end,i)); 

            plot(sigT(j,:,i),'linewidth',4);
%             plot(sigT_run(j,:,i),'r','linewidth',2);
%             plot(sigT_still(j,:,i),'b','linewidth',2);
            plot([prestim prestim],[-1 orientations/2],'--','linewidth',1);% stimulus on time
            title(['cell#' num2str(i) 'ori' num2str(j)]);
            axis([1 seg ymin ymax]);
            axis off
        end
        
%     end
    
end
%% 

% sp=conv(speed,ones(1,prestim)/prestim,'same'); %conv the speed into with a 1s filter
% velocity=interp1(0:.5:info.max_idx, sp(1:info.max_idx*2+1),1:info.max_idx); % downsampling the speed into the frame num
% sigT_run=sigT;
% sigT_still=sigT;
% run=zeros(orientations,rep);% tag speed of each segseg=min(info.frame(3:2:end)-info.frame(1:2:(end-2)))+prestim; % for each segment, prestim, TTL on, TTL off are counted togehter for one seg 
% peak_run=peak;
% peak_still=peak;
%% 


%% 


% figure;title(['cell' num2str(i) 'tuning responses']);hold on
% for i=1:ncell
%     for k=1:orientations
%         subplot(ncell,orientations,(i-1)*orientations+k);
%         plot(sigT(k,:,i));axis tight;
%         %ylim([0 1])
%         peak(k,i)=max(sigT(k,:,i));   %peak is calculated from the 15ms from the stimulus onset towards the end, but each frame is 60ms...
%     end
%     hold off;
% end

    
%%plot a tuning curve and  polar graph based on peak response
margin=.020;
xrow=2;
ycol=ncell;
xpos= margin+((1:xrow)-1)/xrow;
ypos= margin+((1:ycol)-1)/ycol;
xwidth=(1+margin)/xrow-margin;
yheight= (1+margin)/ycol-margin;
ymax=max(peak(:));

y=[peak(:,:); peak(1,:)];
% y_run=[peak_run(:,:); peak_run(1,:)];
% y_still=[peak_still(:,:); peak_still(1,:)];

x=[1:(orientations+1)]';
xp=[0:2*pi/orientations:2*pi]';
xx=[1:.1:(orientations+1)]';

anglesRads = 2*pi*(1:orientations)'/orientations;
for j=1:ncell
gOSI(j) = abs(sum(peak(:,j).*exp(2i*anglesRads)/sum(peak(:,j))));
OSI(j)  = (max(peak(:,j)-min(peak(:,j))))/(max(peak(:,j)+min(peak(:,j))));
end


figure('Name',[fn 'tuning curve &polar graph'],'position',[400 100 1200 1000]);
for i=1:ncell
    subplot('position',[ ypos(i) xpos(1)  yheight xwidth]);
    %try
    %   yy=spline(x,y(:,i),xx);
    %    plot(xx,yy,'-',x,y(:,i),'ok',x,y_run(:,i),'r--',x,y_still(:,i),'b-');
    %catch
    plot(x,y(:,i),'ok-');%,x,y_run(:,i),'r*--',x,y_still(:,i),'b*-'
    %end
    ylabel(['cell' num2str(i)])
    axis([ 1 orientations+1 ymin ymax])
   % axis off
    
    subplot('position',[ ypos(i) xpos(2)  yheight xwidth]);
    P = polar(xp, ymax * ones(orientations+1,1));
    set(P, 'Visible', 'off');
    hold on
%     polar(xp,y_run(:,i),'--r');
%     polar(xp,y_still(:,i),'-b');
    polar(xp,y(:,i),'-k');
    axis off
    
end
    


%%  

%  anglesRads = pi/180*angles;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     osi = abs(sum(FiringRates.*exp(2i*anglesRads)/sum(FiringRates)));
    









