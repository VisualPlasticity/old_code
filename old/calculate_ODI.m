%% Load file
[file,path]=uigetfile('*signals.mat','select contralateral sigF file');
signals1=load(fullfile(path,file))
matrix1 = signals1.matrix;
if isfield(signals1,'window')
    win_sig1 = [signals1.window(1) signals1.window(end)];
else
    win_sig1 = [signals1.win_sig(1) signals1.win_sig(end)];
end
sigF1 = signals1.sigF;
% [~,SI.peakR1,SI.errorR1,~,SI.peakS1,SI.errorS1]= sigFcmp(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
% peak1=fixpeak(SI.peakR1);

[file2,path2]=uigetfile('*signals.mat','select ipsilateral sigF file');
signals2=load(fullfile(path2,file2))
matrix2 = signals2.matrix;
if isfield(signals2,'window')
    win_sig2 = [signals2.window(1) signals2.window(end)];
else
    win_sig2 = [signals2.win_sig(1) signals2.win_sig(end)];
end
sigF2 = signals2.sigF;

[f,p]=uigetfile('*.mat','selected neurons');
load(fullfile(p,f));

%% 
[~,SI.peakR1,SI.errorR1,~,SI.peakS1,SI.errorS1]= sigAcmp(sigF1,win_sig1,matrix1);  % sigR:seg,1,Var,ncell
peak1=fixpeak(SI.peakR1);
[pref1,pref_ori1]=nanmax(peak1(:,1:end-1,:),[],2);

[~,SI.peakR2,SI.errorR2,~,SI.peakS2,SI.errorS2]= sigAcmp(sigF2,win_sig2,matrix2);  % sigR:seg,1,Var,ncell
peak2=fixpeak(SI.peakR2);
[pref2,pref_ori2]=nanmax(peak2(:,1:end-1,:),[],2);

peak0 = (peak1(1,:,pair1)+peak2(1,:,pair2))/2;
[pref0,pref_ori0]=nanmax(peak0(:,1:end-1,:),[],2);

%%
ncell=numel(pair1);
ori_seperate = nan(1,ncell);
ori_combined = nan(1,ncell);
ODI1 = nan(1,ncell);
ODI0 = nan(1,ncell);

%% seperate fitting
nrow=ceil(sqrt(ncell));
ncol=nrow*2;
figure('Position',[0 0 ncol*200 nrow*200])
for i=1:ncell
    subplot(nrow,ncol,i*2-1);hold on
    nth1=pair1(i);
    [temp1,ori1,~]= fit2peakgaussion(pref_ori1(1,:,nth1),peak1(1,1:end-1,nth1));
    text(0,-temp1/10,sprintf('#%d:%.3f',nth1,temp1))
    
    subplot(nrow,ncol,i*2);hold on
    nth2=pair2(i);
    [temp2,ori2,~] = fit2peakgaussion(pref_ori2(1,:,nth2),peak2(1,1:end-1,nth2));
    text(0,-temp2/10,sprintf('#%d:%.3f',nth2,temp2))
    
    if temp1>=temp2
        ori = ori1;
    else
        ori = ori2;
    end
    ODI1(i) = (peak1(1,ori,nth1)-peak2(1,ori,nth2))/(peak1(1,ori,nth1)+peak2(1,ori,nth2));
    ori_seperate(i) = ori ;
    text(14,-temp2/10,sprintf('(%.2f)',ODI1(i)))
end
saveas(gcf,[fullfile(p,f) 'seperate.fig']);
saveas(gcf,[fullfile(p,f) 'seperate.png']);
%% combined fitting
nrow=ceil(sqrt(ncell));
ncol=nrow;
figure('Position',[0 0 ncol*200 nrow*200])
for i=1:ncell
    subplot(nrow,ncol,i);hold on
    [~,ori,~]= fit2peakgaussion(pref_ori0(1,:,i),peak0(1,1:end-1,i));
    text(0,0,sprintf('ori%d',ori))
    
    ori_combined(i) = ori ;
    nth1=pair1(i);
    nth2=pair2(i);
    ODI0(i) = (peak1(1,ori,nth1)-peak2(1,ori,nth2))/(peak1(1,ori,nth1)+peak2(1,ori,nth2));   
    text(14,0,sprintf('(%.2f)',ODI0(i)),'VerticalAlignment','top');
end
saveas(gcf,[fullfile(p,f) 'combined.fig']);
saveas(gcf,[fullfile(p,f) 'combined.png']);
%%
h=figure;
hold on;
seg=.3;
histogram(ODI1,[-1.2:seg:1.2]);
title(sprintf('total cell %d, using preferred orientation', ncell));
ODIfigname= [fullfile(p,f) 'good' num2str(sum(ncell)) '_ODI_seperate'];
saveas(h,[ODIfigname '.fig']);
saveas(h,[ODIfigname '.png']);

h=figure;
hold on;
seg=.3;
histogram(ODI0,[-1.2:seg:1.2]);
title(sprintf('total cell %d, using preferred orientation', ncell));
ODIfigname= [fullfile(p,f) 'good' num2str(sum(ncell)) '_ODI_combined'];
saveas(h,[ODIfigname '.fig']);
saveas(h,[ODIfigname '.png']);
%%
save(fullfile(p,f),'ODI1','seperate','ODI0','ori_combined','-append');

