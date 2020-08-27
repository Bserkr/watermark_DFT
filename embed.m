%基于傅立叶域的数字水印
%注意:水印必须为40*40的二值图像 
%因为40阶的二维arnold置乱周期为30,所以嵌入时置乱8次,提取时置乱22.可以根据自己的需要更改.
%嵌入源码
clc
clear all;

% 保存开始时间
start_time=cputime;
iTimes=8;      %置乱次数
k=1.5;                            % 设置嵌入强度系数
blocksize=8;                    % 块的大小
filter_m=[  1,1,1,1,1,1,1,1;    % 滤波矩阵
            1,1,1,1,1,1,1,1;
            1,1,0,0,0,0,1,1;
            1,1,0,0,0,0,1,1;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;];
        
% 读入载体图像
original_image=imread('lena.jpg');
file_name_g=original_image(:,:,2);%提取绿色分量
cover_object=double(file_name_g)/255;
% 载体图像矩阵的行数与列数
Mc=size(cover_object,1);	        
Nc=size(cover_object,2);	       

% 最大嵌入信息量
max_message=Mc*Nc/(blocksize^2);

% 读入水印图像
mark=imread('mark40.bmp');
mark=rgb2gray(mark);
mark=im2bw(mark);
message=double(mark);
%水印图像矩阵的行数与列数
Mm=size(message,1);	                
Nm=size(message,2);	                

% 检查水印信息是否过大
if  Mm*Nm>max_message   % 载体图像的行*列要大于102400 
   error('水印信息过大')
end

%对水印图像进行Arnold置乱
if Mm~=Nm
    error('水印矩阵必须为方阵');
end
if Mm~=40
    error('必须为40*40大小,或者修改置乱次数');
end
tempImg=message;
for n=1:iTimes % 次数
    for u=1:Mm
        for v=1:Nm
          temp=tempImg(u,v);
          ax=mod(u+v,Mm)+1;
          ay=mod(u+2*v,Nm)+1;
          outImg(ax,ay)=temp;
        end
    end
  tempImg=outImg;
end
message_vector=reshape(outImg,1,Mm*Nm);

% 将cover_object（原图绿色通道矩阵）写入watermarked_image
withmark_image=cover_object;

%置随机数发生器的状态为1100
key=1100;
rand('state',key);

% 产生伪随机序列
pn_sequence_zero=round(2*(rand(1,sum(sum(filter_m)))-0.5));
pn_sequence_one=round(2*(rand(1,sum(sum(filter_m)))-0.5)); 
% 将图像分块
x=1;
y=1;
h=waitbar(0,'嵌入水印，请等待');
for (kk = 1:length(message_vector))
    % 做傅立叶变换
    fft_block=fft2(cover_object(y:y+blocksize-1,x:x+blocksize-1));
    %计算幅值
    abs_block=fftshift(abs(fft_block));
    %计算相位
    angle_block=angle(fft_block);
    % 当message_vector=0且filter_m=1时用伪随机序列pn_sequence_zero叠加abs_block
    % 当message_vector=1且filter_m=1时用伪随机序列pn_sequence_one叠加abs_block
    ll=1;
    if (message_vector(kk)==0)
        for ii=1:blocksize
            for jj=1:blocksize
                if (filter_m(ii,jj)==1)
                    abs_block_o=abs_block(ii,jj);
                    abs_block(ii,jj)=abs_block(ii,jj)*(1+k*pn_sequence_zero(ll));
  abs_block(blocksize-ii+1,blocksize-jj+1)=abs_block(blocksize-ii+1,blocksize-jj+1)+abs_block(ii,jj)-abs_block_o;
                    ll=ll+1;
                end
            end
        end
    else                                     
        for ii=1:blocksize                    
            for jj=1:blocksize
                if (filter_m(ii,jj)==1)
                    abs_block_o=abs_block(ii,jj);
                    abs_block(ii,jj)=abs_block(ii,jj)*(1+k*pn_sequence_one(ll));
  abs_block(blocksize-ii+1,blocksize-jj+1)=abs_block(blocksize-ii+1,blocksize-jj+1)+abs_block(ii,jj)-abs_block_o;
                    ll=ll+1;
                end
            end
        end
    end
    
    % 进行傅立叶逆变换
    abs_block=fftshift(abs_block);
    withmark_image(y:y+blocksize-1,x:x+blocksize-1)=abs(ifft2(abs_block.*exp(i*angle_block)));    
    
    % 移动到下一块
    if (x+blocksize) >= Nc
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
    waitbar(kk/length(message_vector),h);
end
close(h);

% 转换为uint8，并写入
watermarked_image_g=(withmark_image*255);
original_image(:,:,2)=watermarked_image_g;
imwrite(original_image,'withwatermarked_image.bmp','bmp');
%计算运行时间
elapsed_time=cputime-start_time

%计算psnr
PSNR=psnr(cover_object,withmark_image);

% 显示水印，嵌入水印图像与原始图像
figure(1)
subplot(2,2,2)
imshow(message);
title('原始水印');
subplot(2,2,3);
imshow(tempImg);
title('置乱水印');
subplot(2,2,4)
imshow(original_image);
name='嵌入水印图像';
title(strcat(num2str(name),'   k=',num2str(k),'   PSNR=',num2str(PSNR)));
subplot(2,2,1)
imshow(original_image)
title('原始图像');



