%提取
clc
clear all;

% 保存开始时间
start_time=cputime;
iTimes=22;              %置乱次数
blocksize=8;            % 设置块的大小
filter_m=[  1,1,1,1,1,1,1,1;    % 滤波矩阵
            1,1,1,1,1,1,1,1;
            1,1,0,0,0,0,1,1;
            1,1,0,0,0,0,1,1;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;];
% 读入嵌入水印图像
withmark_1=imread('withwatermarked_image.bmp');
withmark=withmark_1(:,:,2);
imshow(withmark_1)
watermarked_image=double(withmark)/255;

% 嵌入水印图像矩阵的行数与列数
Mw=size(watermarked_image,1);	
Nw=size(watermarked_image,2);	        

% 最大可嵌入信息量
max_message=Mw*Nw/(blocksize^2);

% 读入原始水印
mark=imread('mark40.bmp');
mark=rgb2gray(mark);
orig_watermark=double(mark);

% 原始水印矩阵的行数与列数
Mo=size(orig_watermark,1);	
No=size(orig_watermark,2);	

%置随机数发生器的状态为1100
key=1100;
rand('state',key);

% 产生伪随机序列
pn_sequence_zero=round(2*(rand(1,sum(sum(filter_m)))-0.5));
pn_sequence_one=round(2*(rand(1,sum(sum(filter_m)))-0.5));              
% 将图像分块
x=1;
y=1;
h=waitbar(0,'提取水印，请等待');
for (kk = 1:max_message)

    % 傅立叶变换
    fft_block_w=fft2(watermarked_image(y:y+blocksize-1,x:x+blocksize-1));
    abs_block_w=abs(fftshift(fft_block_w));
    ll=1;
    for ii=1:blocksize
        for jj=1:blocksize
            if (filter_m(ii,jj)==1)
                sequence(ll)=abs_block_w(ii,jj); 
                ll=ll+1;
            end
        end
    end
   
    % 计算sequence与pn_sequence_zero和pn_sequence_one的相关系数
   correlation_zero(kk)=corr2(pn_sequence_zero,sequence);
   correlation_one(kk)=corr2(pn_sequence_one,sequence);               
    
    % 移动到下一块
    if (x+blocksize) >= Nw
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
    waitbar(kk/max_message,h);
end
close(h);

% 如果correlation_zero>correlation_one，那么message_vector=0，反之为1.
for (kk=1:Mo*No)
    if correlation_zero(kk)>correlation_one(kk)
        message_vector(kk)=0;
    else
        message_vector(kk)=1;
    end
end

% Arnold置乱
tempImg=reshape(message_vector(1:Mo*No),Mo,No);
message_arnold=tempImg;
for n=1:iTimes % 次数
    for u=1:Mo
        for v=1:No
          temp=tempImg(u,v);
          ax=mod(u+v,Mo)+1;
          ay=mod(u+2*v,No)+1;
          outImg(ax,ay)=temp;
        end
    end
  tempImg=outImg;
end
message=outImg;
message=uint8(message);
imwrite(outImg,'marked.jpg','jpg')

% 计算运行时间
elapsed_time=cputime-start_time,
%计算NC（归一化相关系数）
NC=nc(message,orig_watermark)

% 显示提取水印与原始水印
figure(3)
subplot(1,2,1);
imshow(message,[]);
name='提取水印';
title(strcat(num2str(name),'NC=',num2str(NC)));
subplot(1,2,2)
imshow(orig_watermark,[])
title('原始水印');