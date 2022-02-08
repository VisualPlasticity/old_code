function ball = sbxballspeed(fn)


load(fn,'time','ball');

dT=diff(time);
speed=abs(ball(1:end-1))/192*2./dT'; %size of the imaging area 192 pixal->2cm

fname=fn(1:(findstr(fn,'_ball')-1));
save([fname '.ball'],'dT','speed');

figure,plot(Interleave(time,time(2:end-1)),Interleave(speed,speed));title (fn);

[value,pos]=max(speed);
nframe=25;
pos=max(floor(pos/2),nframe+1); %find max speed points and display -nframe-+nframe for varification
colormap gray;
load(fn,'data');
figure;
for i=(pos-nframe):(pos+nframe)
imshow(data(:,:,i));
text(0,0,sprintf('frame %8d',i),'color',[1 1 1],'fontsize',10);
drawnow;
end