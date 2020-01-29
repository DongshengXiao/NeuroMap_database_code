function [SM] = OIA_STA(I,Jc,SF,refper,win)
% function [SM] = OIA_STA(I,Jc,SF,refper,win)
% input:
% I = [X,Y,F] imaging
% J = timestamp (s)
% SF = sampling freq (Hz)
% win = temporal window (how many sec before and after trig), eg: [-5 5]
% refper = refractory period [min max] (sec), eg: [3 4]  if negative, apply both way !!!!!!!
% output: % SM = [X,Y,F,T] of pull

if refper(2)<0
    refper = abs(refper);
    tw = 1;
else
    tw = 0;
end

Jc = Jc(Jc>0); % no zero and no <0
good = [];
if tw == 0
    for i = 2:length(Jc)
        d = Jc(i)-Jc(i-1);
        if ((d>refper(1))&(d<refper(2)))
            good = [good i];
        end
    end
else
    for i = 2:length(Jc)-1
        d = Jc(i)-Jc(i-1);
        d2 = Jc(i+1)-Jc(i);
        if (((d>refper(1))&(d<refper(2)))&((d2>refper(1))&(d2<refper(2))))
            good = [good i];
        end
    end
end
Jc = Jc(good);

win = round(win.*SF);
winsize = win(2)-win(1)+1;

SM = NaN.*ones(size(I,1),size(I,2),winsize,length(Jc)-2,'single'); 
for i = 2:length(Jc)-1
    loca = round(Jc(i)*SF);
    disp(['pull ' num2str(i) '/' num2str(length(Jc)-1) ' (t=' num2str(Jc(i)) 's)'])
    try
        S = single(I(:,:,loca+win(1):loca+win(2)));
         
        if win(1)<0
            M = mean(S(:,:,1:abs(win(1))),3);
        else
            M = mean(S(:,:,1:size(S,3)),3);
        end
        for f = 1:size(S,3)
            S(:,:,f) = 100*(S(:,:,f) - M)./M;
        end
        SM(:,:,:,i) = S;  
    catch
        disp(' > out of range')
    end
end

smm = squeeze(mean(mean(mean(SM,1),2),3));
SM = SM(:,:,:,~isnan(smm));
