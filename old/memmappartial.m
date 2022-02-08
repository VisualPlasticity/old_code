function fname=memmap(nam,magnification,picsiz);
%% load file
%path_to_package = 'C:\Users\kraken\Documents\MATLAB\ca_source_extraction\utilities';   % path to the folder that contains the package
%addpath(genpath(path_to_package));
%close all;

starttime=datestr(now);
data = matfile([nam '_memmap.mat']);
if nargin<2
    try
        magnification=data.magnification;
    catch
        magnification=2;
    end
end


%% Set parameters
try
    sizY=data.sizY;
catch
    sizY = size(data,'Y');                  % size of data matrix
end
if nargin<3
    if numel(strfind(nam,'_'))<3
        picsiz=[1 sizY(1) 1 796] + [5 -5 5 -5];
    else
        picsiz=[25 sizY(1) 10 sizY(2)-10];%two planes bidirectional
    end
end
pcsX = picsiz(2)-picsiz(1) ;
pcsY = picsiz(4)-picsiz(3);
overlap = floor([4,6]*magnification);                        % amount of overlap in each dimension (optional, default: [4,4])
square=80;
patch_size(1) = ceil(mod(pcsX,square)/floor(pcsX/square))+square ;
patch_size(2) = ceil(mod(pcsY,square)/floor(pcsY/square))+square ;
% size of each patch along each dimension (optional, default: [128,120])
patches = construct_partial_patches(picsiz,patch_size+overlap,overlap);
K=ceil(20/magnification^2/150^2*square^2);
tau = floor([3,5]*magnification);          % std of gaussian kernel (half width/height of neuron)
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.6;                                  % merging threshold
tsub=20;
options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'ssub',magnification,...
    'tsub',tsub,...
    'save_memory',1,...
    'max_size', ceil(6*magnification),...
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );
%     'fudge_factor',0.98,...                     % bias correction for AR coefficients

%display(sprintf('CNMFparms done %d min',round(toc/60)));

%% Run on patches
tic;

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

display(sprintf('run patches in %d min',round(toc/60)));

%% DONnotorder , just plot
contour_threshold = 0.95;                 % amount of energy used for each component to construct contour plot
K_m = size(C,1);
%[A_or,C_or,S_or,P] = order_ROIs(A,C,S,P); % order components
%     [sig,S_df] = extract_DF_F(data,[A_or,b],[C_or;f],K_m+1); % extract DF/F values (sig) and spiking info
    [sig,S_df] = extract_DF_F(data,[A,b],[C;f],K_m+1); % extract DF/F values (sig) and spiking info

% try
%     V=data.V;
% catch
%     m=info.aligned.m;
%     V=zeros(sizY(1),sizY(2));
%     fac=sqrt(sizY(end));
%     for i=1:sizY(end)
%         V=V+((double(data.Y(:,:,i)-m))/fac).^2;
%     end
% end
V=reshape(P.sn,sizY(1),sizY(2));

%    [Coor,json_file] = plot_contours(A_or,V,contour_threshold,1);
    [Coor,json_file] = plot_contours(A,V,contour_threshold,1);
%% save signals into different files
resultname=sprintf('%s_%dcell_tau%d_%d_tsub%d',nam,K_m,tau(1),tau(2),tsub);
fname=strrep(resultname,'\\mps-zfs\data\jsun','\\mps-pc53\');
sig=sig';
save([fname '.signals'],'sig','S_df');
savejson('jmesh',json_file,[fname '.jmesh']);        % optional save json file with component coordinates (requires matlab json library)
saveas(gcf,[fname '.png'],'png');
try 
    saveas(gcf,[fname '.fig'],'fig');
catch
end
try
    T=0;
    for n=1:numel(data.eachsize)
        sig_chunk=sig(T+1:T+data.eachsize(1,n),:);
        save([fname '_' num2str(n) '.signals'],'sig_chunk','S_df');
        T=T+data.eachsize(1,n);
    end
catch
end
endtime=datestr(now);
display(fname)
display(starttime);
display(endtime);
% sendmail('j.suninchina@gmail.com','Finished project', ...
%          [fname 'finished']);