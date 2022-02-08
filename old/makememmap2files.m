function newname=makememmapnfiles;
%% define the name and save to MPS-ZFS directly
path=pwd;
zfs_path=strrep(path,'c:','\\mps-zfs\data\jsun');

files = dir('*.sbx');
newname=[files(1).name(1:end-5) 'x_memmap.mat'];
newname=fullfile(zfs_path,newname);

if exist(newname,'file')
    disp([newname ' memmap file already existing']);return;
elseif ~exist(zfs_path,'dir')
    mkdir(zfs_path);
end
data = matfile(newname,'Writable',true);
disp(sprintf(' try saving memmap file: %s',newname));

%% save each piece
global info;            % this contains the information about the structure of the image
%calculate the total size of the image
T0=0;
for i=1:numel(files)
    fn=files(i).name(1:end-4);
    sbxread(fn,1,1);% read one frame to read the header of the image sequence
    T0=T0+info.max_idx+1;
end

T=0;
for i=1:numel(files)
    fn=files(i).name(1:end-4);
    sbxread(fn,1,1);% read one frame to read the header of the image sequence
    m(:,:,i)=info.aligned.m;
    eachsize(i)=info.max_idx;
    
    V=zeros(info.sz);
    fac=sqrt(info.max_idx);
    tic
    if i==1
        data.m=zeros(info.sz);
        data.Y=zeros([info.sz T0],'uint16');
        data.V=zeros(info.sz);
        u=0;v=0;
    else
        [u v] = fftalign(m(:,:,i),m(:,:,1));
        figure('Name',['aligned file#' num2str(i)]);
        imshowpair(circshift(squeeze(m(:,:,i)),[u v]),squeeze(m(:,:,1)))
    end
    for j=0:info.max_idx
        z =sbxread(fn,j,1);
        z= squeeze(z);
        img = circshift(z,info.aligned.T(j+1,:)+[u v]); % align the image
        data.Y(:,:,T+j+1) = img;
        V=V+((double(img-m(:,:,i)))/fac).^2;
        if mod(j,500)==0
            fprintf('File %d Frame %d/%d for %.2f seconds\n ',i,j,info.max_idx,toc);
        end
    end
    save([fn '.align'],'V','-append');
    ratio=T/(T+eachsize(i));
    data.V=data.V*ratio+circshift(V,[u v])*(1-ratio);
    data.m=data.m*ratio+circshift(double(m(:,:,i)),[u v])*(1-ratio)
    T=T+eachsize(i);
end
%assert(T==info.max_idx,'the total image size is not matched!')
data.eachsize=eachsize;
data.sizY =[info.sz T0];
data.Yr = reshape(data.Y,prod(info.sz),T0);
data.nY = min(reshape(data.Yr,prod(info.sz)*T0,1));
if mod(info.config.magnification*2,5) == 0
    data.magnification = info.config.magnification*2/5;
else
    data.magnification = info.config.magnification;
end

