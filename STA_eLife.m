
%% 1. load ePhys and image
clear all
fileePhys = 'C:\Users\user\Dropbox\nmpaper\Figure 1-source data 1\2015-5-6-1_cell_1.txt'; % spike time file
fileImg = 'C:\Users\user\Dropbox\nmpaper\Figure 1-source data 1\2015-5-6-1.tif'; % imaging file
 
sp = dlmread(fileePhys); % open ePhys file
 
I = OIA_wfi_load_simple(fileImg); % open Imaging
I = imresize(I,1/2); % binning 2x2

SF = size(I,3)./(sp(end)-sp(1)); % calculate sampling frequency *Hz)
disp(['Sampling Frequency: ' num2str(SF) ' Hz'])

sp = sp-sp(1); sp(1)=[]; % spike train

%% 2. average with pulling trigger
clc
win = [-3 3]; % temporal window (how many sec before and after trig)
refper = [0 Inf]; % refractory period [min max] (sec)

spr = sort(max(sp)*rand(1,length(sp))); % random spike vector

[SM] = OIA_STA(I,sp,SF,refper,win); % generate STA matrix from spikes
[SMr] = OIA_STA(I,spr,SF,refper,win); % ... and random spikes

SMM = mean(SM,4); % average over trials
SMMr = mean(SMr,4);

 
%% 3. check signal
k=1; % roi size
 
map = max(SMM,[],3) + min(SMM,[],3); % maximum and minimum T projection map
range = max(abs([max(map(:)) min(map(:))])); range = [-range range]; % zero centered range
x=k+1; y=k+1;
while 1
    subplot(121), imshow(map,range), colormap jet; colorbar
    hold on, plot(x,y,'wo'); plot(x,y,'xk'); hold off
    
    s = squeeze(mean(mean(SM(y-k:y+k,x-k:x+k,:,:),1),2)); % signal at location x,y for each trial
    s_s = squeeze(mean(mean(SMr(y-k:y+k,x-k:x+k,:,:),1),2));
    
    ms = mean(s,2); % calculating mean and SEM
    ss = std(s,[],2)./sqrt(size(s,2));   
    mss = mean(s_s,2);
    sss = std(s_s,[],2)./sqrt(size(s_s,2));
    
    tt = ([1:length(ss)]'-round(abs(win(1)).*(SF))-1)/(SF); % time vector
    subplot(122), 

    plot(tt,ms,'k','LineWidth',2); hold on, plot(tt,ms-ss,'k'); plot(tt,ms+ss,'k'); 
    plot(tt,mss,'LineWidth',2,'Color',[.7 .7 .7]); plot(tt,mss-sss,'Color',[.7 .7 .7]); plot(tt,mss+sss,'Color',[.7 .7 .7]);
    plot([min(tt) max(tt)],[0 0],'k')
    range2 = 1*(max(ms)-min(ms));
    hold off 
    grid on, xlabel('time (s)'); ylabel('DF/F (%)')
    title([num2str(size(SM,4)) ' spikes'])    
    
    [x,y]=myginput;
end

