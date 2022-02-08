function [peak,idx,error]=calER(sigF) 
% sigF:seg,rep,Var,ncell
% peak 1*Var*ncell
%error:1,Var*ncell
rep=size(sigF,2);
Var=size(sigF,3);
ncell=size(sigF,4);
sigA=squeeze(mean(sigF,2));
[peak,idx]=max(sigA(window,:,:),[],1);  % peak 1*Var*ncell
temp=permute(sigF,[2 1 3 4 ]);  %temp:rep,seg,Var,ncell
error=std(temp(:,idx),0,1)/sqrt(numel(rep)); %error:1,Var*ncell
error=reshape(error,1,Var,ncell);
end