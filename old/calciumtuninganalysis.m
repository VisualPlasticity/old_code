% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. preprocess:
%   align the file sbxalignxdir
%   manuelly find the cells sbxsegmentpoly
%   match the stmilus
% 	pull the signals sbxpullsignalsdir  
% OR sbxprocess
%   align the file sbxalignxdir
%   find the cells sbxsegmentpoly
% 	pull the signals sbxpullsignalsdir  
%(2) plot signal sbxplotlines
% 3. calculate the ball motion sbxballmotion &eyemotion
% 4. match speed and stimulus data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
% load calcium .mat data
[f,p]=uigetfile('*.mat','load calcium imaging data');
f=fullfile(p,f);
load(f)
fname=f(1:(findstr(f,'.mat')-1));
sbxread(fname,1,1);
global info;
movement=preprocess(fname);
load([fname '.signals'],'-mat')
if movement
try
    load([fname '.ball'],'-mat')
catch
    load([fname '_ball'],'-mat')
end
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%extra and analyze data
ncell=size(sig,2);
for k=1:ncell  
    base=prctile(sig(info.frame(1):end,k),10); % use the ones with stimulus presented as the
    sigF(:,k)=sig(:,k)/base-1;
end
%sig(:,:)=sig(:,:)./sig(:,1)*(1:ncell);

CAframeHz =info.resfreq/info.recordsPerBuffer;
CAlineHz = info.config.lines;
prestim=ceil(CAframeHz*2); %use 2s baseline

orientations=numel(unique(info.stimtype));% calculate types of orientations
rep=length(info.stimtype)/orientations; %calculate repetitions 
sp=conv(speed,ones(1,prestim)/prestim,'same'); %conv the speed into with a 1s filter
velocity=interp1(0:.5:info.max_idx, sp(1:info.max_idx*2+1),1:info.max_idx); % downsampling the speed into the frame num
seg=min(info.frame(3:2:end)-info.frame(1:2:(end-2))); % info.frame contains both TTL ON and TTL off signal, for each segment, TTL on+ TTL off are counted togehter for average 
stimON=info.frame(2)-info.frame(1); % StimON duration

frame_on=zeros(orientations,rep); %tag frame_on of each seg
run=zeros(orientations,rep);% tag speed of each seg
run_rep=zeros(orientations);% try find if there is too much bias? 

sig_seg=zeros(orientations,rep,prestim+seg,ncell);
sigT=zeros(orientations,prestim+seg,ncell); % build sigT based on signals from sig and organized according to segment
sigT_run=zeros(orientations,prestim+seg,ncell);
sigT_still=zeros(orientations,prestim+seg,ncell);

peak=zeros(orientations,ncell); %  
peak_run=peak;
peak_still=peak;

figure('Name','data with stim &ball movement', 'position',[200 100 2400 1000]);
stim=zeros(1,info.max_idx); %just for plotting purpose,assign stimtype to all the ON frames 
for i=1:(size(info.frame,1)/2)
    stim(info.frame(i*2-1):info.frame(i*2))=info.stimtype(i);
end
ymax=max(max(sigF(:)),stim);
ymin=min(sigF(:));
axis_x=1:info.max_idx;
col=1;
row=ncell+movement
[xpos,ypos,xwidth,yheight]=figurepara(col,row,axis_x(1),axis_x(end),ymin,ymax);
for k=1:ncell
    subplot('position',[xpos ypos(k) xwidth yheight]);hold on; 
    plot(axis_x,sigF(:,k),'b','linewidth',2);
    plot(axis_x,stim,'k','linewidth',1);
    axis([axis_x(1) axis_x(end) ymin ymax]);
    text(1,ymax/2,['cell' num2str(k)],'HorizontalAlignment','center');
    axis off
end
if movement
    subplot('position',[xpos ypos(k+1) xwidth yheight]);hold on;    
    plot(axis_x,velocity,'b',axis_x,stim,'k','linewidth',1);
    axis tight;
    plot( [0;0],[0;10],'-k','Linewidth',2);
    text(-1,5, '10cm/s')
    axis off
end


for j=1:orientations
    ori=find(info.stimtype==j);
    on=ori*2-1;
    frame_on(j,:)=info.frame(on);
    for k=1:rep
        run(j,k)=mean(velocity(frame_on(j,k):frame_on(j,k)+stimON-1));%calucate running speed during stimON
        baseline=min(sigF(frame_on(j,k)-prestim:frame_on(j,k),:),[],1); %calculate baseline using 2s prestim
        sig_seg(j,k,seg+prestim,:)=sigF(frame_on(j,k)-prestim:(frame_on(j,k)+seg-1),:)-baseline*(1:seg)'; %prestim +seg 
    
    end
end
sigT=mean(sig_seg,2);
sigT_run=mean(sig_seg(run>=.5,:,:),2);
sigT_stil=mean(sig_seg(run<.5,:,:),2);
peak=max(sigT(:,prestim+1:end,:),[],2);   %peak is calculated from the 15ms from the stimulus onset towards the end, but each frame is 60ms...
peak_run(:,:)=max(sigT_run(:,prestim+1:end,:),[],2);
peak_still(:,:)=max(sigT_still(:,prestim+1:end,:),[],2);

figure('Name','reorganized','position',([ 200 100 1600 1200]));
xrow=orientations;
ycol=ncell;
axis_x=1:seg+prestim;
ymax=max(sigF(:));
ymin=min(sigF(:));
[xpos,ypos,xwidth,yheight]=figurepara(col,row,axis_x(1),axis_x(end),ymin,ymax);
for i=1:ncell
        for j=1:orientations
            subplot('position',[xpos(j) ypos(i) xwidth yheight]);hold on;
            plot(sig_seg(run(j,:)>=.5,:,i),'r--');
            plot(sig_seg(run(j,:)<.5,:,i),'b--');
            plot(sigT(j,:,i),'linewidth',4);
            plot(sigT_run(j,:,i),'r','linewidth',2);
            plot(sigT_still(j,:,i),'b','linewidth',2);
            plot([prestim prestim],[-1 orientations/2],'--','linewidth',1);% stimulus on time
            title(['cell#' num2str(i) 'ori' num2str(j)]);
            axis([1 seg+prestim ymin ymax]);
            axis off
        end 
end

   
%%plot a tuning curve and  polar graph based on peak response
figure('Name','tuning curve &polar graph','position',[400 100 1200 1000]);
xrow=ncell;
ycol=2;
ymax=max(peak(:));
ymin=min(peak(:));

x=[1:(orientations+1)]';
xp=[0:2*pi/orientations:2*pi]';
xx=[1:.1:(orientations+1)]';
y=[peak(:,:); peak(1,:)];
yy=spline(x,y(:,i),xx);
y_run=[peak_run(:,:); peak_run(1,:)];
y_still=[peak_still(:,:); peak_still(1,:)];

[xpos,ypos,xwidth,yheight]=figurepara(col,row,x(1),x(end),ymin,ymax);

for i=1:ncell
    subplot('position',[  xpos(i) ypos(1) xwidth yheight ]);
    % plot(xx,yy,'-',x,y(:,i),'ok',x,y_run(:,i),'r--',x,y_still(:,i),'b-');
    plot(x,y(:,i),'ok',x,y_run(:,i),'r*--',x,y_still(:,i),'b*-');
    title(['cell' num2str(i)])
    axis([ 1 orientations+1 ymin ymax])
    axis off
    
    subplot('position',[  xpos(i) ypos(2) xwidth yheight ]);   
    P = polar(xp, ymax * ones(orientations+1,1));
    set(P, 'Visible', 'off');
    hold on
    polar(xp,y_run(:,i),'--r');
    polar(xp,y_still(:,i),'-b');
    polar(xp,y(:,i),'-k');
    axis off
end
    

% 
% x=x(1:nSteps1);
% osi=((sum(finalvalue.*sin(2*x)))^2+(sum(finalvalue.*cos(2*x)))^2)^0.5/sum(finalvalue)  %
%  anglesRads = pi/180*angles;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     osi = abs(sum(FiringRates.*exp(2i*anglesRads)/sum(FiringRates)));
    










