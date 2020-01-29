function [roi,MC,MC1,sortie] = OIA_roi(varargin)
% fonction [roi,MC,MC1,s] = OIA_roi(argument)
% Matthieu Vanni - 15 mars 2006
% argument obligatoire :
% 'map' : matrice sous la forme [x y z]
%         si une seule carte z = 1
%         sinon les analyses sont moyennees sur tout les z
% arguments optionnels
% 'roi' : image binaire de la roi
% 'resize' : rajuste l'image en %
% 'FS' : parametre de l'analyse DANS L'ORDRE !      
%             [lag roigauss fftr smooth surface largeur]
%         lag = course de l'autocorrelation
%         roigauss = multiplication de l'autocorrelogramme par une
%                    gaussienne avant la TF (cutoff en % de la ROI)
%         fftr = binning de la OITF
%         smooth = lissage de l'autocorrelogramme (en pixel)
%         surface = normalisation par la surface reelle de correlation
%                  (1 = actif)
%         largeur : largeur de l'image en mm pour permettre la mesure en mm
%         par defaut :
%             lag = 1000 (quand lag depasse la taille de la ROI,
%                         lag prend automatiquement la taille de la
%                         ROI)
%             roigauss = 0
%             fftr = 10
%             smooth = 0
%             surface = 1
%             VECTEUR : [20 0 10 0 1]
%         3 formes d'entree :
%             [lag roigauss fftr smooth surface largeur]
%             [lag]
%             []
% exemples :
% OIA_roi('map',I,'resize',0.5,'FS',[20 .5 10 3 1 13])
% OIA_roi('map',I,'FS',[20])
% OIA_roi('map',I,'FS',[])
%
% * * * SPECIAL MARTIN * * * SPECIAL MARTIN * * * 
% [roi,A,B] = OIA_roi('martin',M,R,Z)
% [roi,A,B] = OIA_roi('martin2',M,R,Z)
% avec M = matrice 3D de quantif
%      R = matrice 2D de "clicage" ou binaire
%      Z = coordonnees en z
% si Z = 1 -> 1 histogramme
% si Z = 2 -> 2 histogrammes + courbe m/s + scatter
% dans ce cas avec [A,B] vecteur du scatter
% et s est le signal moyen pour chaque condition
% si Z > 2 -> courbe m/s
% Histogramme : medianne en vert, moyenne en rouge
% si 'martin' = ROI par clicage
% si 'martin2' = ROI pre-determinee dans la matrice binaire R
% * * * SPECIAL MARTIN * * * SPECIAL MARTIN * * * 
%
% --------------------
% implementer comme OIA_clic et avoir la possibilite d'avoir 3 conditions
% (ex avant/pendant/apres)
MC = []; % TRUC A LA CON A ENLEVER !
MC1 = []; % TRUC A LA CON A ENLEVER !
roi = [];
entree = varargin;
indice = OIA_reco_entree(entree);
param = OIA_stru_entree(entree,indice);
param = OIA_resize(param);
param = OIA_setroi(param);
roi = param.BW;
if sum(sum(param.quantif)) ~= 0
    [param,sortie] = OIA_quantif(param);
    MC = param.S;
    MC1 = param.S1;
else
    [X,Y,Z]=size(param.I);
    if Z == 2
        param = OIA_crosscorr(param);
    else
        param = OIA_autocorr(param);
    end;
    try, param = OIA_fft(param); end
    try, param = OIA_poolCorr(param); end
    try, param = OIA_davis(param); end
    try, OIA_afficheFS(param); end
    try, MC = param.MC; end
end;

function [param,m] = OIA_quantif(param);
[X,Y,Z]=size(param.quantif);
[A,B]=sort(param.BW(:),1);
for i = 1:Z
    K = param.quantif(:,:,i);
    K = K(:);
    S = K(B(min(find(A == 1)):length(A)));
    % ou sont les nans !
    wnan = find(isnan(S) == 1);
    if wnan ~= 0
        debut = 1;
        for j = 1:length(wnan)
            S2 = S(debut : wnan(j) - 1);
            debut = wnan(j) + 1;
        end;
        S = S2;
    end;
    m(i) = mean(S);
    s(i) = std(double(S));
    if i == 1, S1 = S; end;
end
if Z == 1
    texte = ['m=',num2str(mean(S)),', md=',num2str(median(S)),', s=',num2str(std(S))];
    [param.S,param.S1] = hist(S,100);
    hist(S,100),title(texte)
    line([mean(S) mean(S)],[0 max(hist(S,100))],'Color',[1 0 0])
    line([median(S) median(S)],[0 max(hist(S,100))],'Color',[0 1 0])
elseif Z == 2
    texte = ['m=',num2str(mean(S1)),', md=',num2str(median(S1)),', s=',num2str(std(S1))];
    subplot(221), hist(S1,100), title(texte)
    line([mean(S1) mean(S1)],[0 max(hist(S1,100))],'Color',[1 0 0])
    line([median(S1) median(S1)],[0 max(hist(S1,100))],'Color',[0 1 0])
    texte = ['m=',num2str(mean(S)),', md=',num2str(median(S)),', s=',num2str(std(S))];
    subplot(222), hist(S,100), title(texte)    
    line([mean(S) mean(S)],[0 max(hist(S,100))],'Color',[1 0 0])
    line([median(S) median(S)],[0 max(hist(S,100))],'Color',[0 1 0])
     
    subplot(223), 
    param.S = S;
    param.S1 = S1;
%     H1 = S - min(S);
%     H1 = 100 * H1 / max(H1);
%     H2 = S1 - min(S1);
%     H2 = 100 * H2 / max(H2); 
%     OIA_correlation(H1,H2); colormap jet, title('NORMALISE EN %')
    
    scatter(S1,S,'.')
    
    subplot(224), errorbar(param.abs,m,s);
else
    if param.spectre == 1
        plot(param.abs,m)
    elseif param.spectre == 0
        try, errorbar(param.abs,m,s); end
    end;
end;

function param = OIA_davis(param);
K = param.K;
K = K > mean(K(:));
imagesc(K);
indice = 1;
for i = 1:2*param.LX+1
    for j = 1:2*param.LY+1
        if K(i,j) == 1
            D(indice) = sqrt((((param.LY+1)-j)^2) + (((param.LX+1)-i)^2));
            indice = indice + 1;
        end;
    end;
end;
[param.davis.H]=hist(D,20);
param.davis.X = max(D)*(1:20)/20;

function param = OIA_resize(param);
if param.resize ~= 1
    disp(['resize = ',num2str(param.resize)])
    if param.quantif ~= 0
        [X Y Z] = size(param.quantif);
        for i = 1:Z
            K(:,:,i) = imresize(param.quantif(:,:,i),param.resize);
        end;
        param.clic = imresize(param.clic,param.resize);
        param.quantif = K;       
    else
        for i = 1:param.nmap
            K(:,:,i) = imresize(param.I(:,:,i),param.resize);
        end;
        param.I = K;
    end;
end;

function param = OIA_poolCorr(param);
for i = 1:360
    for j = 1:param.LX
        A = 2*pi*i/360;
        dX = round(param.LX+1+j*cos(A));
        dY = round(param.LX+1+j*sin(A));
        MC(i,j) = param.K(dX,dY);
    end;
end;
param.MC = MC(1:180,:) + MC(181:360,:); 
for i = 1:360
    for j = 1:param.LX*param.fftr
        A = 2*pi*i/360;
        dX = round(param.fftr*param.LX+1+j*cos(A));
        dY = round(param.fftr*param.LX+1+j*sin(A));
        MF(i,j) = param.TF(dX,dY);
    end;
end;
param.MF = MF(1:180,:) + MF(181:360,:); 
  
function OIA_afficheFS(param)
[X,largeur]=size(param.I);
mmpp = param.L /largeur; % mm par pixel
R = param.I(:,:,1) + param.BW*max(max(param.I(:,:,1)));
subplot(241), imagesc(R)
title('ROI')
subplot(242), 
if mmpp == 0
    imagesc(param.K)
    xlabel('distance (pixel)'), ylabel('distance (pixel)'), title('autocorrelation')
else
    limite = length((param.K-1)/2) * mmpp / 2;
    imagesc([-limite limite],[-limite limite],param.K)
    xlabel('distance (mm)'), ylabel('distance (mm)'), title('autocorrelation')
end;
XX = (-param.fftr*param.LX:+param.fftr*param.LX)/(param.fftr*2*param.LX+1);
subplot(243), imagesc(XX,XX,param.TF)
xlabel('frequence (cycle/pixel)'), ylabel('frequence (cycle/pixel)'), title('TF(autocorrelation)')
subplot(245), 
if mmpp == 0
    imagesc(param.MC)
    xlabel('distance (pixel)'), ylabel('axe (degre)'), title('profil de correlation')
else
    [X,Y] = size(param.MC);
    imagesc([0:mmpp*Y],[],param.MC);
    xlabel('distance (mm)'), ylabel('axe (degre)'), title('profil de correlation')
end;
subplot(246), 
if mmpp == 0
    s = sum(param.MC);
    s = 100 * s / max(s);
    plot(s)
    xlabel('distance (pixel)'), ylabel('correlation %'), title('somme des correlations')
    axis([0 param.LX min(s) max(s)])
else
    [X,Y] = size(param.MC);
    dy = 1/(Y-1);
    y = [0:dy:1]*mmpp*Y;
    s = sum(param.MC);
    s = 100 * s / max(s);
    plot(y,s)
    xlabel('distance (mm)'), ylabel('correlation %'), title('somme des correlations')
    axis([0 max(y) min(s) max(s)])
end;



XX = (1:+param.fftr*param.LX)/(param.fftr*2*param.LX+1);
YY = [0:180];
subplot(247), imagesc(XX,YY,param.MF)
xlabel('frequence (cycle/pixel)'), ylabel('axe (degre)'),
title('profil de la TF(correlation)')
subplot(248), plot(XX,sum(param.MF))
xlabel('frequence (cycle/pixel)'), ylabel('amplitude de la TF'),
title('somme de la TF(correlation)')
axis([0 max(XX) 0 max(sum(param.MF))])
subplot(244), plot(param.davis.X,param.davis.H)
xlabel('distance (pixel)'), ylabel('comptage'),
title('distance eucl. (Meth. Davis 1990)')
axis([0 max(param.davis.X) 0 max(param.davis.H)])

function param = OIA_fft(param)
if param.surface == 1
    K = param.C ./ param.S;
else
    K = param.C;
end;
if param.roigauss ~= 0
    K = fspecial('gaussian',size(K),length(K)*param.roigauss) .* K;
end;
if param.smooth ~= 0
    masque = fspecial('gaussian',param.smooth*4,param.smooth);
    K = conv2(K,masque,'same');
end;
param.TF = abs(fft2(K,2*param.LX*param.fftr+1,2*param.LY*param.fftr+1));
param.TF = fftshift(param.TF);
param.K = K;

function param = OIA_crosscorr(param);
disp('crosscorrelation');
I = param.I(:,:,1) .* param.BW;
I = I(min(param.y):max(param.y),min(param.x):max(param.x));
J = param.I(:,:,2) .* param.BW;
J = J(min(param.y):max(param.y),min(param.x):max(param.x));
if max(param.x) - min(param.x) < param.LX
    param.LX = max(param.x) - min(param.x);
end;
if max(param.y) - min(param.y) < param.LY
    param.LY = max(param.y) - min(param.y);
end;
if param.LY > param.LX
    param.LY = param.LX;
else
    param.LX = param.LY;
end;
[X Y] = size(I);
REF = zeros(X+2*param.LX,Y+2*param.LY);
REF(param.LX+1:param.LX+X,param.LY+1:param.LY+Y) = J;
for i = -param.LX:param.LX
    disp(['correlation / lag = ', num2str(i)]);
    for j = -param.LY:param.LY
        LAG = zeros(X+2*param.LX,Y+2*param.LY);
        LAG(param.LX+1+i:param.LX+X+i,param.LY+1+j:param.LY+Y+j) = I;     
        C(i+param.LX+1,j+param.LY+1) = sum(sum(REF.*LAG));
        S(i+param.LX+1,j+param.LY+1) = (X-abs(i))*(Y-abs(j));
    end;
end;
param.C = mean(C,3);
param.S = mean(S,3);

function param = OIA_autocorr(param);
for m = 1:param.nmap
    disp(['------Correlation sur la carte ', num2str(m)]);
    I = param.I(:,:,m) .* param.BW;
    I = I(min(param.y):max(param.y),min(param.x):max(param.x));
    if max(param.x) - min(param.x) < param.LX
        param.LX = max(param.x) - min(param.x);
    end;
    if max(param.y) - min(param.y) < param.LY
        param.LY = max(param.y) - min(param.y);
    end;
    if param.LY > param.LX
        param.LY = param.LX;
    else
        param.LX = param.LY;
    end;
    [X Y] = size(I);
    REF = zeros(X+2*param.LX,Y+2*param.LY);
    REF(param.LX+1:param.LX+X,param.LY+1:param.LY+Y) = I;
    for i = -param.LX:param.LX
        disp(['correlation / lag = ', num2str(i)]);
        for j = -param.LY:param.LY
            LAG = zeros(X+2*param.LX,Y+2*param.LY);
            LAG(param.LX+1+i:param.LX+X+i,param.LY+1+j:param.LY+Y+j) = I;     
            C(i+param.LX+1,j+param.LY+1,m) = sum(sum(REF.*LAG));
            S(i+param.LX+1,j+param.LY+1,m) = (X-abs(i))*(Y-abs(j));
        end;
    end;
end;
try, param.C = mean(C,3); end
try, param.S = mean(S,3); end


function param = OIA_setroi(param);
if param.setroi == 1
    if param.quantif == 0
        [Y X Z] = size(param.I);
        for i = 1:Z
            param.I(:,:,i) = param.I(:,:,i) - min(min(param.I(:,:,i)));
            param.I(:,:,i) = 2*(param.I(:,:,i) ./max(max(param.I(:,:,i))))-1;
        end;
        R = param.I(:,:,1);
    else
        R = param.clic;
        [Y X Z] = size(R);
    end;
    imagesc(R); colormap gray
    i = 1; xi = 2; yi = 2; x = 0;
    while (xi > 1 & yi > 1 & xi < X & yi < Y)  
        [xi,yi] = ginput(1);
        if (xi > 1 & yi > 1 & xi < X & yi < Y)
            x(i) = round(xi);
            y(i) = round(yi);
            if length(x) >= 2
                line([x(i-1) x(i)],[y(i-1) y(i)],'Color',[1 0 0],'LineWidth',1);
            end;
            i = i + 1;
        end;
    end;
    if length(x) <3
        x = [1 X X 1];
        y = [1 1 Y Y];
    else
        x(i) = x(1);
        y(i) = y(1);
    end;
    param.BW = roipoly(R,x,y);
    param.x = x;
    param.y = y;
else
    disp('>>>>> Matrice Binaire en entree')
    disp('>>>>> Pas de definition manuelle de la ROI')
    param.BW = param.clic;
end;

function indice = OIA_reco_entree(entree)
indice.resize = 0;
indice.fs = 0;
indice.map = 0;
indice.martin = 0;
indice.martin2 = 0;
indice.spectre = 0;
i = 1;
while i <= length(entree)
    resu = cell2struct(entree(i),{'donnee'},2);
    if ischar(resu.donnee) == 1
        switch resu.donnee
            case 'resize', indice.resize = i; 
            case 'FS', indice.fs = i;
            case 'map', indice.map = i;
            case 'martin', indice.martin = i;
            case 'martin2', indice.martin2 = i;
            case 'spectre', indice.spectre = i;
        end;
    end;
    i = i + 1;
end;

function param = OIA_stru_entree(entree,indice)
param.setroi = 1;
if indice.map ~= 0
    resu = cell2struct(entree(indice.map + 1),{'donnee'},2);
    param.I = double(resu.donnee);
    [X Y param.nmap] = size(param.I);
else
    param.I = 0;
    param.nmap = 0;
end;
if indice.resize ~= 0
    resu = cell2struct(entree(indice.resize + 1),{'donnee'},2);
    param.resize = resu.donnee;
else
    param.resize = 1;
end;
if indice.fs ~= 0
    resu = cell2struct(entree(indice.fs + 1),{'donnee'},2);
    data = resu.donnee;
    if length(data) == 1
        param.LX = data; % facultatif, mettre au max et ca fittera avec la ROI
        param.LY = data;
        param.roigauss = 0; % en pourcentage de la ROI de corr
        param.fftr = 1; 
        param.smooth = 0; % en pourcentage de la ROI de corr
        param.surface = 1;
        param.L = 0;
    elseif length(data) > 1
        param.LX = data(1);
        param.LY = data(1);
        param.roigauss = data(2);
        param.fftr = data(3); 
        param.smooth = data(4);
        param.surface = data(5);
        param.L = data(6);
    elseif length(data) == 0
        param.LX = 10000;
        param.LY = 10000;
        param.roigauss = 0;
        param.fftr = 1; 
        param.smooth = 0;
        param.surface = 1;
        param.L = 0;
    end;
end;
param.quantif = 0;
if indice.martin ~= 0
    resu = cell2struct(entree(indice.martin + 1),{'donnee'},2);
    param.quantif = resu.donnee;
    resu = cell2struct(entree(indice.martin + 2),{'donnee'},2);
    param.clic = resu.donnee;
    resu = cell2struct(entree(indice.martin + 3),{'donnee'},2);
    param.abs = resu.donnee;
    param.setroi = 1;
    param.spectre = 0;
end;
if indice.martin2 ~= 0
    resu = cell2struct(entree(indice.martin2 + 1),{'donnee'},2);
    param.quantif = resu.donnee;
    resu = cell2struct(entree(indice.martin2 + 2),{'donnee'},2);
    param.clic = resu.donnee;
    resu = cell2struct(entree(indice.martin2 + 3),{'donnee'},2);
    param.abs = resu.donnee;
    param.setroi = 0;
    param.spectre = 0;
end;
if indice.spectre ~= 0
    resu = cell2struct(entree(indice.spectre + 1),{'donnee'},2);
    param.quantif = resu.donnee;
    resu = cell2struct(entree(indice.spectre + 2),{'donnee'},2);
    param.clic = resu.donnee;
    resu = cell2struct(entree(indice.spectre + 3),{'donnee'},2);
    param.abs = resu.donnee;
    param.setroi = 1;
    param.spectre = 1;
end;
param.S = 0;
param.S1 = 0;
