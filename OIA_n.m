function I = OIA_n(I)
% function I = OIA_n(I);
% I = I - min(I(:));
% I = I / max(I(:));
for i = 1:size(I,3)
    for j = 1:size(I,4)
        for k = 1:size(I,5)
            Ii = I(:,:,i,j,k);
            Ii = Ii - min(Ii(:));
            I(:,:,i,j,k) = Ii / max(Ii(:));
        end
    end
end
