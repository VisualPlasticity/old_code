function sbxalimemdir(s)

if(nargin>=1) % cell with filenames to be aligned
    d=dir([s '.sbx']);
else
    d = dir('*.sbx');
end

% Align all *.sbx files in the list

for(i=1:length(d))
    try
        fn = strtok(d(i).name,'.');
        if(exist([fn '_memmap.mat'])==0) 
            sbxread(fn,1,1);        % read one frame to read the header of the image sequence
            global info Y;                % this contains the information about the structure of the image
            tic
            [m,T] = sbxalimem(fn,0:info.max_idx-1);   %
            sizY = size(Y);
            Yr = reshape(Y,prod(sizY(1:end-1)),[]);
            nY = min(Yr(:));
            Yr = Yr - nY;
            save([fn '.align'],'m','T');
            save([fn '_memmap.mat'],'Yr','Y','nY','sizY','-v7.3');
            display(sprintf('Done %s: Aligned %d images and saved memmapin %d min',fn,info.max_idx,round(toc/60)));
            clear Y Yr info;
%             h=figure('Name',fn);
%             imagesc(m);
%             savefig(h,[fn num2str(info.max_idx)]);
        else
            sprintf('File %s is already aligned and made memmap file',fn)
        end
    catch
        sprintf('Could not make align&makememmap %s',fn)
    end
end