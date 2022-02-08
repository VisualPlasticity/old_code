function ca_analysis(f);
% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;global info;

%% load calcium signal and plot
load(f,'-mat') %f should be the .signals

if exist([f(1:end-8) '.fig'])
    openfig([f(1:end-8) '.fig'])
end

%% load stimulus data and plot
fname=f(1:findstr(f,'_memmap')-1)
movement=preprocess(fname);
CAframeHz =info.resfreq/info.recordsPerBuffer;
CAlineHz = info.config.lines;

if size(sig,1)<size(sig,2)
    sig=sig';  
end

sig=sig(:,1:end-1);  %%%% only plot HALF of the cells auto-selected

stim=zeros(1,info.max_idx); %just for plotting purpose,assign stimtype to all the ON frames
for i=1:(numel(info.frame)/2)
    stim(info.frame(i*2-1):info.frame(i*2))=info.stimtype(i);
end

sigplt(sig,stim/max(stim)); %

%% set paramenter for stimulus
ncell=size(sig,2);

prestim=ceil(CAframeHz*1); %use 1s baseline
stimON=info.frame(2)-info.frame(1); % StimON duration
seg=prestim+min(info.frame(3:2:end)-info.frame(1:2:(end-2))); % segment=prestim+info.frame(TTL ON+TTL off )

Var=prod(info.steps);% calculate types of stimulus, orientations, contrast, etc
rep=floor(min(2*numel(info.stimtype),numel(info.frame))/2/Var); %calculate repetitions

sigT=zeros(seg*rep,Var,ncell); % build sigT based on signals from sig and organized according to segment
frame_on=zeros(Var,rep); %tag frame_on of each seg
peak=zeros(Var,ncell); %
%% reorganize stimulus according to stimulus type
for j=1:Var
    ori=find(info.stimtype==j);
    on=ori*2;%-1;
    frame_on(j,:)=info.frame(on(1:rep));%-prestim-5; % the frame marking the prestim timepoint size 1*rep
    temp=ones(seg,1)*frame_on(j,:)+ (0:seg-1)'*ones(1,rep); % array seg*rep
    sigT(:,j,:)=sig(temp,:);
end
sigF=reshape(sigT,seg,rep,Var,ncell); % sigF:seg,rep,Var,ncell
sigA=mean(sigF,2); % sigA seg*Var*ncell
sigA=squeeze(sigA);
peak=max(sigA(prestim:end,:,:),[],1);
peak=squeeze(peak);

sigFplt(sigF,prestim);  % sigF:seg,rep,Var,ncell
%%
%% syn with ball/eye motion
if movement
    load([fname '.ball'],'-mat');
    load([fname '.eye'],'-mat');

    sp=conv(speed,ones(1,prestim)/prestim,'same'); %conv the speed into with a 1s filter
    velocity=interp1(0:1/2:info.max_idx, sp(1:info.max_idx*2+1),1:info.max_idx); % downsampling the speed into the frame num
    figure;hold on;plot(stim);plot(velocity);
    
    run=[];%frame_on=zeros(Var,rep);
    figure;
    for i=1:stimON+prestim
        run=cat(3,run,velocity(frame_on+i-1));
        imagesc(run(:,:,i));
        drawnow;
    end
    figure;subplot(2,1,1);imagesc(mean(run,3));colorbar; subplot(2,1,2);imagesc(max(run,[],3));colorbar;
     
    % sigF:seg,rep,Var,ncell
    [run_trial(:,2),run_trial(:,1)]=find(squeeze(mean(run,3))>0.1);run_trial
    [still_trial(:,1),still_trial(:,2)]=find(squeeze(mean(run,3))<0.02);
    sigR=mean(sigF(:,run_trial(:,1),run_trail(:,2),:),2);
   % sigS=
end
    %%
   
    %             peak(j,i)=max(sigT(j,prestim+1:end,i));   %peak is calculated from the 15ms from the stimulus onset towards the end, but each frame is 60ms...
    %             peak_run(j,i)=max(sigT_run(j,prestim+1:end,i));
    %             peak_still(j,i)=max(sigT_still(j,prestim+1:end,i));
    
    %             plot(sigT_run(j,:,i),'r','linewidth',2);
    %             plot(sigT_still(j,:,i),'b','linewidth',2);
    %             peak_run(j,i)=max(sigT_run(j,prestim+1:end,i));
    %             peak_still(j,i)=max(sigT_still(j,prestim+1:end,i));
    %             plot(sigT_run(j,:,i),'r','linewidth',2);
    %             plot(sigT_still(j,:,i),'b','linewidth',2);
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
    
    
    %% 
    %%plot a tuning curve and  polar graph based on peak response
    margin=.020;
    xcol=2;
    yrow=ncell;
    ymax=max(peak(:));
    [xpos,ypos,xwidth,yheight]=figurepara(xcol,yrow,xmin,xmax,ymin,ymax)
    
    y=[peak(:,:); peak(1,:)];
    % y_run=[peak_run(:,:); peak_run(1,:)];
    % y_still=[peak_still(:,:); peak_still(1,:)];
    
    Ori=Var;  %
    x=[1:(Ori+1)]';
    xp=[0:2*pi/Ori:2*pi]';
    xx=[1:.1:(Ori+1)]';
    
    anglesRads = 2*pi*(1:Ori)'/Ori;
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
        axis([ 1 Ori+1 ymin ymax])
        % axis off
        
        subplot('position',[ ypos(i) xpos(2)  yheight xwidth]);
        P = polar(xp, ymax * ones(Ori+1,1));
        set(P, 'Visible', 'off');
        hold on
        %     polar(xp,y_run(:,i),'--r');
        %     polar(xp,y_still(:,i),'-b');
        polar(xp,y(:,i),'-k');
        axis off
        
    end
    %% 
    
    
    
    
    
    
    
    
    
    
    
