clc;
clear;
close  all
I = imread('color.bmp');    %��ȡ����ͼ��
I = rgb2gray(I);         %ת��Ϊ�Ҷ�ͼ
W = imread('mark128.bmp');    %��ȡˮӡͼ��
W = rgb2gray(W); 
figure('Name','����ͼ��')
imshow(I);
title('����ͼ��')
figure('Name','ˮӡͼ��')
imshow(W);
title('ˮӡͼ��')
ntimes=23;    %��Կ1�� Arnold���Ҵ���
rngseed=59433;     %��Կ2�����������
flag=1;      %�Ƿ���ʾ�м�ͼ��
[Iw,psnr]=setdwtwatermark(I,W,ntimes,rngseed,flag);    %ˮӡǶ��
[Wg,nc]=getdwtwatermark(Iw,W,ntimes,rngseed,flag);      %ˮӡ��ȡ
%%  �����Ҫ���й������ԣ��������⼸�д���ȡ��ע�͡�
% close all    
% action={'filter','resize','crop','noise','rotate'};      %ˮӡ����
% for i = 1:numel(action)
%     dwtwatermarkattack(action{i},Iw,W,ntimes,rngseed);
% end

