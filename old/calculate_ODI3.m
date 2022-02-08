% For one figure

clear all;close all;clc
%
%% load data and fig
[figname1,path1a]=uigetfile('.fig','find the fig');

[~,ipsieyepath]=uigetfile('*.mat','select trace for ipsieyepath');
[~,contraeyepath]=uigetfile('*.mat','select trace for contraeyepath');
% caanalysisplot(ipsieyepath)
% caanalysisplot(contraeyepath)

% [f1,p1]=uigetfile('*.signals','select .signals for ipsieyepath');
% [f2,p2]=uigetfile('*.signals','select .signals for contraeyepath');
% ipsieyepath=caanalysis(fullfile(p1,f1));
% contraeyepath=caanalysis(fullfile(p2,f2));
%%
eye1=load(fullfile(ipsieyepath,'peakSI.mat'));
if isfield(eye1,'window')
    win_sig1 = [eye1.window(1) eye1.window(end)];
else
    win_sig1 = [eye1.win_sig(1) eye1.win_sig(end)];
end
if ~isfield(eye1,'matrix')
    eye1.matrix=eye1.run.matrix;
end
if isempty(eye1.matrix)
    eye1.matrix = logical(ones(size(eye1.sigF,2),size(eye1.sigF,3)));
end
if ~isfield(eye1.SI,'errorS')
    eye1.SI.peakR=eye1.SI.peak;
    eye1.SI.errorR = eye1.SI.error;
    eye1.SI.peakS=eye1.SI.peak;
    eye1.SI.errorS = eye1.SI.error;
end
eye1.sigF= permute(eye1.sigF,[4 3 2 1]);
eye1.finalvalueR = permute(eye1.SI.peakR,[3 2 1]);
eye1.stdofeachvalueR = permute(eye1.SI.errorR,[3 2 1]);
eye1.finalvalueS = permute(eye1.SI.peakS,[3 2 1]);
eye1.stdofeachvalueS = permute(eye1.SI.errorS,[3 2 1]);

eye2=load(fullfile(contraeyepath,'peakSI.mat'));
if isfield(eye2,'window')
    win_sig2 = [eye2.window(1) eye2.window(end)];
else
    win_sig2 = [eye2.win_sig(1) eye2.win_sig(end)];
end
if ~isfield(eye2,'matrix')
    eye2.matrix=eye2.run.matrix;
end
if ~isfield(eye2.SI,'errorS')
    eye2.SI.peakR=eye2.SI.peak;
    eye2.SI.errorR = eye2.SI.error;
    eye2.SI.peakS=eye2.SI.peak;
    eye2.SI.errorS = eye2.SI.error;
end
eye2.sigF= permute(eye2.sigF,[4 3 2 1]);
eye2.finalvalueR = permute(eye2.SI.peakR,[3 2 1]);
eye2.stdofeachvalueR = permute(eye2.SI.errorR,[3 2 1]);
eye2.finalvalueS = permute(eye2.SI.peakS,[3 2 1]);
eye2.stdofeachvalueS = permute(eye2.SI.errorS,[3 2 1]);
[~,pref_ori1]=nanmax(eye1.finalvalueR(:,1:end-1,:),[],2);
[~,pref_ori2]=nanmax(eye2.finalvalueR(:,1:end-1,:),[],2);

%% Fig info
fig1=openfig(fullfile(path1a,figname1));
fig1img =getimage(fig1);
fig1img = fig1img/max(fig1img(:));
kC1= fig1.Children.Children;
ncell=(numel(kC1)-1)/2;
SIZ = size(fig1img);
F1Contour =zeros(prod(SIZ),ncell);
u=0;v=0;
for i=1:ncell
    nth = ncell+1-i;
    F1cell(i).Position = [kC1(nth).Position(1)+v+5,kC1(nth).Position(2)+u+5];
    F1cell(i).String = kC1(nth).String;
    F1Contour(:,i) = reshape(circshift(kC1(nth+ncell).ZData,[u,v]),prod(SIZ),1);
end

%% Define parameters
pref_ori = nan(1,ncell);
ORI = nan(1,ncell);
ORI0= nan(1,ncell);
ORI1= zeros(1,ncell);
ODI = nan(1,ncell);
ODI1 = nan(1,ncell);
ODI0= nan(1,ncell);
ODI_bi = nan(1,ncell);
ODI1_bi = nan(1,ncell);
ODI0_bi= nan(1,ncell);
nSteps = size(eye1.finalvalueR,2);
framestocapture=min (size(eye1.sigF,4), size(eye2.sigF,4));

%% Select cells based on responses and image
% maskpic1 = zeros(SIZ);
[f,p]=uigetfile('*.mat','Load region from previous file?');
if f~= 0
    load(fullfile(p,f),'maskpic1');
else
    maskpic1=roipoly;
end
pair0 = find(reshape(maskpic1,1,prod(SIZ))*F1Contour >0);
%% Select cells based on responses

h0=figure;
imshow(fig1img)
hold on;
%235:numel(pair0)
for k=1:numel(pair0)
    j= pair0(k);
    nth1= j;
    nth2= j;
%     h0;
%     contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1])
%     text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[ 1 1 1],'FontSize',8)
%     drawnow;

    peak1 = eye1.finalvalueR(nth1,:);
    peak2 = eye2.finalvalueR(nth2,:);
    [a1,ori1]=max(peak1);
    [a2,ori2]=max(peak2);
    [a,pref_ori]=max(cat(2,peak1,peak2));
    pref_ori = mod(pref_ori,nSteps);
    [odi,odi_bi] = ODIcalculation(pref_ori,peak2,peak1);
    
    h1=figure('Position',[100 200 1500 400],'Name',['alltraces_cell#' num2str(j)]);
    subplot(3,3,1);    hold on;
    errorbar(1:nSteps,eye1.finalvalueS(nth1,:),eye1.stdofeachvalueS(nth1,:),'Linewidth',1,'Color',[0 0 0 0.3]);
    errorbar(1:nSteps,eye1.finalvalueR(nth1,:),eye1.stdofeachvalueR(nth1,:),'Linewidth',2,'Color',[0 1 0 0.3]);
    
    %     xlim([1 nSteps]);
    axis tight;
    title(['ipsi eye, ori=' num2str(ori1)] );
    
    subplot(3,3,2);    hold on;
    errorbar(1:nSteps,eye2.finalvalueS(nth2,:),eye2.stdofeachvalueS(nth2,:),'Linewidth',1,'Color',[0 0 0 0.3]);
    errorbar(1:nSteps,eye2.finalvalueR(nth2,:),eye2.stdofeachvalueR(nth2,:),'Linewidth',2,'Color',[1 0 0 0.3]);
    %     xlim([1 nSteps]);
    axis tight;
    title(['contra eye,ori=' num2str(ori2)] )
    
    subplot(3,3,3);    hold on;
%     plot(1:nSteps,eye1.finalvalueS(nth1,:),'Linewidth',2,'Color',[0 1 0 0.3]);
    plot(1:nSteps,eye1.finalvalueR(nth1,:),'Linewidth',2,'Color',[0 1 0 0.7]);
%     plot(1:nSteps,eye2.finalvalueS(nth2,:),'Linewidth',2,'Color',[1 0 0 0.3]);
    plot(1:nSteps,eye2.finalvalueR(nth2,:),'Linewidth',2,'Color',[1 0 0 0.7]);
    xlim([1 nSteps]);
%     legend('ispi-R','ispi-S','contra-R','contra-S' );
    legend('ispi-R','contra-R');
    legend('boxoff','Location','top') ;
    title(sprintf('a=%.2f,ori=%d,odi=%.2f,odi_b_i=%.2f',a,pref_ori,odi,odi_bi));
    
    temp_1st=squeeze(eye1.sigF(nth1,:,:,:));
    temp_2nd=squeeze(eye2.sigF(nth2,:,:,:));
    ymax=prctile([temp_1st(:);temp_2nd(:)],99);
    ymin=0;
     buffer=(ymax-ymin)*.05;
    for ori=1:nSteps
        traces1 = squeeze(temp_1st(ori,eye1.matrix(:,ori),:))';
        traces1s = squeeze(temp_1st(ori,~eye1.matrix(:,ori),:))';
%         Fr1=pwrplt(nanmean(traces1,2));

        h1;b=subplot(3,nSteps,nSteps+ori); hold on;
        plot(traces1s,'Color',[0 0 0 0.3]);
        plot(traces1,'Color',[0 1 0 0.3]);
        plot(nanmean(traces1s,2),'Linewidth',1,'Color',[0 0 0 0.7]);
        plot(nanmean(traces1,2),'Linewidth',2,'Color',[0 1 0 0.7]);
        plot([win_sig1(1) win_sig1(1)],[ymin ymax],'Linewidth',1,'Color',[0 0 0 0.3]);% stimulus on time
        plot([win_sig1(end) win_sig1(end)],[ymin ymax],'Linewidth',1,'Color',[0 0 0 0.3]);% stimulus off time
        title(sprintf('%.3f',eye1.finalvalueR(nth1,ori)));
%         title(sprintf('%.2f',Fr1))
        b.XLabel.String=sprintf('%d',ori);
        axis([1 framestocapture ymin-buffer ymax+buffer]);
        
        traces2= squeeze(temp_2nd(ori,eye2.matrix(:,ori),:))';
        traces2s= squeeze(temp_2nd(ori,~eye2.matrix(:,ori),:))';
%         Fr2=pwrplt(nanmean(traces2s,2));
        
        h1;c=subplot(3,nSteps,2*nSteps+ori); hold on;
        plot(traces2s,'Color',[0 0 0 0.3]);
        plot(traces2,'Color',[1 0 0 0.3]);
        plot(nanmean(traces2s,2),'Linewidth',1,'Color',[0 0 0 0.7]);
        plot(nanmean(traces2,2),'Linewidth',2,'Color',[1 0 0 0.7]);
        plot([win_sig2(1) win_sig2(1)],[ymin ymax],'Linewidth',1,'Color',[0 0 0 0.3]);
        plot([win_sig2(end) win_sig2(end)],[ymin ymax],'Linewidth',1,'Color',[0 0 0 0.3]);
        axis([1 framestocapture ymin-buffer ymax+buffer]);
        title(sprintf('%.3f',eye2.finalvalueR(nth2,ori)));
% title(sprintf('%.2f',Fr2))

        if ori==1
            b.YLabel.String=sprintf('ipsieye,cell# %d',nth1);
            c.YLabel.String=sprintf('contraeye,cell# %d',nth2);
             b.XTickLabel='';

        else
            b.YTickLabel='';
            b.XTickLabel='';
            c.YTickLabel='';
            c.XTickLabel='';
        end
    end
    % %     m=inputdlg('Good cell& best ori=','Selected cell?',1,{num2str(pref_ori)})
    % %     if ~isempty(m)
    % %         ORI(j) =str2num(m{1});
    ORI1(j) = pref_ori;
    PEAK1(j)= a;
    ODI1(j) = odi;
    ODI1_bi(j) = odi_bi;
    saveas(h1,[h1.Name '.fig'])
    saveas(h1,[h1.Name '.png'])
% %         else
    % %          ORI(j) =0;
    % %         ORI1(j) =0;
    % %     end
%     close(h1)
%     h0;
%     if isnan(odi)
%         text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[ 1 1 1],'FontSize',8)
%         contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1])
%     else
%         text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[odi+(odi<0) 1-odi-(odi<0) 0],'FontSize',8)
%         contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[odi+(odi<0) 1-odi-(odi<0) 0])
%         %      text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[ odi>0 odi<0 0],'FontSize',8)
%         %      contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[odi>0 odi<0 0])
%     end
%         drawnow;

end

% fign=strtok(figname1,'t');
% saveas(h0,[fign 'selected.fig']);
% saveas(h0,[fign 'selected.png']);

%%

% ORI1(isnan(ORI1))=0;
% ORI1(ORI1==17)=0;
% pair1=find(ORI1~=0);
% for k=1:numel(pair1)
%     j=pair1(k);
%         nth1= j;
%         nth2= j;
% %         [ODI(j),ODI_bi(j)]= ODIcalculation(ORI(j),eye2.finalvalueR(nth2,:),eye1.finalvalueR(nth1,:));
%         [ODI1(j),ODI1_bi(j)]= ODIcalculation(ORI1(j),eye2.finalvalueR(nth2,:),eye1.finalvalueR(nth1,:));
% %         [ODI0(j),ODI0_bi(j)]= ODIcalculation(ORI0(j),eye2.finalvalueR(nth2,:),eye1.finalvalueR(nth1,:));
% end

% %%
h2=figure;
imshow(fig1img)
hold on;
%contra eye,red ;ipsi eye,green
for m=1:numel(pair1)
    j= pair1(m);
    if isnan(odi)
        text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[ 1 1 1],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1])
    else
        text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[ODI1(j)+(ODI1(j)<0) 1-ODI1(j)-(ODI1(j)<0) 0],'FontSize',8)
        contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[ODI1(j)+(ODI1(j)<0) 1-ODI1(j)-(ODI1(j)<0) 0])
        %      text(F1cell(j).Position(1),F1cell(j).Position(2),F1cell(j).String,'Color',[ odi>0 odi<0 0],'FontSize',8)
        %      contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,j),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[odi>0 odi<0 0])
    end
end

% fign=strtok(figname1,'t');
% saveas(h2,[fign 'selected.fig']);
% saveas(h2,[fign 'selected.png']);
%%
pair1=find(ORI1~=0);
ODI1_a=nanmean(ODI1(pair1));
ODI1_b=nanmean(ODI1_bi(pair1));
figtitle=sprintf('Total %d cell,ODI=%.2f', numel(pair1),ODI1_a);
h=figure;
seg=.3;
subplot(1,2,1);hold on;title(sprintf('ODI-direction =%.2f',ODI1_b));
histogram(ODI1(pair1),[-1.2:seg:1.2]);
subplot(1,2,2);hold on;title(figtitle);
histogram(ODI1_bi(pair1),[-1.2:seg:1.2]);

saveas(h,[fign 'ODI.fig']);
saveas(h,[fign 'ODI.png']);
%%
save([fign 'ODI.mat'],'ODI1','ODI1_bi','PEAK1','ODI1_a','pair0','pair1','maskpic1');

