close all;
%clear all; clc
load magnetization.txt
figure(1);
for i=1:250000
imagesc(magnetization(1+(i-1)*200:i*200,:));axis equal;axis square;colorbar
%title(strcat("T=",num2str(i*0.01)))
title(i);
%saveas (1, strcat("T",num2str(i*0.01),".png"));
drawnow;pause(.3)
end

