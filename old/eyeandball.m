%function eyeandball(fn);
%clearvars -except fn;
clear all;

fn='bal_818_034';
preprocess(fn);
showtime=200;
sbxread(fn,1,1);
global info
%% read in stimulus ON, 
stim=zeros(1,info.max_idx); %just for plotting purpose,assign stimtype to all the ON frames
for i=1:(numel(info.frame)/2)
    stim(info.frame(i*2-1):info.frame(i*2))=info.stimtype(i);
end

%% read in stimulus, eye, ball, mousetrack

CAframeHz =info.resfreq/info.recordsPerBuffer;
sec=ceil(CAframeHz);

tsamp=info.mousetrack.tsamp;
MT=info.mousetrack.vsmooth;


load([fn '.ball'],'ball','dT','-mat');
sp=abs(ball(2:end))/192*2.5./dT'; %size of the imaging area,192pixel=2.5cm
sp=conv(sp,ones(1,sec)/sec,'same'); %conv the speed into with a 1s filter
bsamp=cumsum(dT)*CAframeHz;
speed=interp1(bsamp,sp,1:numel(bsamp)/2);

%load([fn '.eye'],'Area','-mat');
%A=conv(Area,ones(1,sec)/sec,'same'); 



[~,pos]=max(speed);
pos=max(pos,showtime*2);
%   
% load([fn '_eye.mat'],'data');
% data = squeeze(data); % the raw images...
% data_eye=data(:,:,pos-showtime:pos+showtime);
% clear data;

% load([fn '_ball.mat'],'data');
% data = squeeze(data);
% data_ball=data(:,:,(pos-showtime)*2:2:(pos+showtime)*2);
% clear data;




%% do the plotting
%clf;

%figure;
subplot(2,1,1);title('mouse1')
hold on;
%plot(tsamp,speed);
plot(1:numel(bsamp)/2,speed)
plot(tsamp,MT);  %+info.frame(1)
%plot((A-prctile(A,5))/(max(A(500:2000))-prctile(A,5)));
legend(['ball' num2str(sum(speed)/CAframeHz)],['mousetrack' num2str(sum(MT)/CAframeHz)]);%eye
%xlim([0 2500])
drawnow;

% colormap gray;
% for i=1:showtime*2+1
% subplot(2,2,3);imshow(data_eye(:,:,i));
% subplot(2,2,4);imshow(data_ball(:,:,i));
% drawnow;
% end
%% 




