clc;
clear;
close  all
I = imread('color.bmp');    %读取载体图像
I = rgb2gray(I);         %转换为灰度图
W = imread('mark128.bmp');    %读取水印图像
W = rgb2gray(W); 
figure('Name','载体图像')
imshow(I);
title('载体图像')
figure('Name','水印图像')
imshow(W);
title('水印图像')
ntimes=23;    %密钥1， Arnold置乱次数
rngseed=59433;     %密钥2，随机数种子
flag=1;      %是否显示中间图像
[Iw,psnr]=setdwtwatermark(I,W,ntimes,rngseed,flag);    %水印嵌入
[Wg,nc]=getdwtwatermark(Iw,W,ntimes,rngseed,flag);      %水印提取
%%  如果需要进行攻击测试，则将下面这几行代码取消注释。
% close all    
% action={'filter','resize','crop','noise','rotate'};      %水印攻击
% for i = 1:numel(action)
%     dwtwatermarkattack(action{i},Iw,W,ntimes,rngseed);
% end

