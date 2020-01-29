function [x,y,bn]=myginput
% function [x,y,bn]=myginput
[x,y,bn]=ginput(1);
x = round(x);
y = round(y);
