function [M,roi,SF,bx,by] = OIA_wfi_load_simple(file)
%function [M,roi,SF,bx,by] = OIA_wfi_load_simple(file)


info = imfinfo(file);
tampon = imread(file,'Index',1);
F = length(info);
M = zeros(size(tampon,1),size(tampon,2),F,'uint16');
M(:,:,1) = tampon(:,:,1);
tic
wait_bar = waitbar(0,'open file');
ind = 0;
for i = 2:F
    if ind == 0, waitbar(i/F, wait_bar); end;
    ind = ind + 1; if ind == 100, ind = 0; end;
    tampon = imread(file,'Index',i,'Info',info);
    M(:,:,i) = tampon(:,:,1);
end;
close(wait_bar);
temps = num2str(round(10*toc)/10);
disp([file ' open in ' num2str(temps) 's'])


fileroi = [file(1:length(file)-4) '_roi.mat'];
resu = dir(fileroi);
if length(resu) == 0
    
    s = squeeze(sum(sum(M,1),2));
    plot(s); ylim([0 max(s)])
    [x,y]=ginput(1);
    s = bwlabel(s>y);
    M = M(:,:,(s==1));
    
    roid = OIA_n(single(M(:,:,round(size(M,3)/2)))).^.5;
    
    imshow(roid), colormap gray, title('BREGMA ????')
    [bx,by]=myginput;

    
    title('>>>> ROI LEFT'); roi1 = OIA_roi('martin',roid,roid,1);
    title('>>>> ROI RIGHT'); roi2 = OIA_roi('martin',roid,roid,1);
    roi = (roi1 + roi2)>0;
    close all
    
    SF = input('what is the sampling frequency?');
    rot = input('need a rotation?: how much degree?');
    M = imrotate(M,rot,'crop'); roi = imrotate(roi,rot,'crop');
    roid = single(M(:,:,1));

    imagesc(roid), colormap gray, 
    title('where is anterior?'); [x,y] = myginput;
    if y>size(M,1)/2, 
        fud = 1; 
        disp('flip upside down'); 
        M = M([size(M,1):-1:1],:,:); 
        roi = roi([size(roi,1):-1:1],:); 
        roid = single(M(:,:,1));
        imagesc(roid), colormap gray
    else, fud = 0; end;
    
    title('where is right side?'); [x,y] = myginput;
    if x<size(M,2)/2, 
        flr = 1; 
        disp('flip left right'); 
        M = M(:,[size(M,2):-1:1],:); 
        roi = roi(:,[size(roi,2):-1:1]); 
    else, flr = 0; end;
    save(fileroi,'roi','fud','flr','rot','SF','s','bx','by')
else
    load(fileroi)
    M = M(:,:,(s==1));
    if rot~=0, disp('rotation'); M = imrotate(M,rot,'crop'); end  
    if fud==1, disp('flip upside down'); M = M([size(M,1):-1:1],:,:);  end;
    if flr==1, disp('flip left right'); M = M(:,[size(M,2):-1:1],:); end;
end
close all




