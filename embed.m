%���ڸ���Ҷ�������ˮӡ
%ע��:ˮӡ����Ϊ40*40�Ķ�ֵͼ�� 
%��Ϊ40�׵Ķ�άarnold��������Ϊ30,����Ƕ��ʱ����8��,��ȡʱ����22.���Ը����Լ�����Ҫ����.
%Ƕ��Դ��
clc
clear all;

% ���濪ʼʱ��
start_time=cputime;
iTimes=8;      %���Ҵ���
k=1.5;                            % ����Ƕ��ǿ��ϵ��
blocksize=8;                    % ��Ĵ�С
filter_m=[  1,1,1,1,1,1,1,1;    % �˲�����
            1,1,1,1,1,1,1,1;
            1,1,0,0,0,0,1,1;
            1,1,0,0,0,0,1,1;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;];
        
% ��������ͼ��
original_image=imread('lena.jpg');
file_name_g=original_image(:,:,2);%��ȡ��ɫ����
cover_object=double(file_name_g)/255;
% ����ͼ����������������
Mc=size(cover_object,1);	        
Nc=size(cover_object,2);	       

% ���Ƕ����Ϣ��
max_message=Mc*Nc/(blocksize^2);

% ����ˮӡͼ��
mark=imread('mark40.bmp');
mark=rgb2gray(mark);
mark=im2bw(mark);
message=double(mark);
%ˮӡͼ����������������
Mm=size(message,1);	                
Nm=size(message,2);	                

% ���ˮӡ��Ϣ�Ƿ����
if  Mm*Nm>max_message   % ����ͼ�����*��Ҫ����102400 
   error('ˮӡ��Ϣ����')
end

%��ˮӡͼ�����Arnold����
if Mm~=Nm
    error('ˮӡ�������Ϊ����');
end
if Mm~=40
    error('����Ϊ40*40��С,�����޸����Ҵ���');
end
tempImg=message;
for n=1:iTimes % ����
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

% ��cover_object��ԭͼ��ɫͨ������д��watermarked_image
withmark_image=cover_object;

%���������������״̬Ϊ1100
key=1100;
rand('state',key);

% ����α�������
pn_sequence_zero=round(2*(rand(1,sum(sum(filter_m)))-0.5));
pn_sequence_one=round(2*(rand(1,sum(sum(filter_m)))-0.5)); 
% ��ͼ��ֿ�
x=1;
y=1;
h=waitbar(0,'Ƕ��ˮӡ����ȴ�');
for (kk = 1:length(message_vector))
    % ������Ҷ�任
    fft_block=fft2(cover_object(y:y+blocksize-1,x:x+blocksize-1));
    %�����ֵ
    abs_block=fftshift(abs(fft_block));
    %������λ
    angle_block=angle(fft_block);
    % ��message_vector=0��filter_m=1ʱ��α�������pn_sequence_zero����abs_block
    % ��message_vector=1��filter_m=1ʱ��α�������pn_sequence_one����abs_block
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
    
    % ���и���Ҷ��任
    abs_block=fftshift(abs_block);
    withmark_image(y:y+blocksize-1,x:x+blocksize-1)=abs(ifft2(abs_block.*exp(i*angle_block)));    
    
    % �ƶ�����һ��
    if (x+blocksize) >= Nc
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
    waitbar(kk/length(message_vector),h);
end
close(h);

% ת��Ϊuint8����д��
watermarked_image_g=(withmark_image*255);
original_image(:,:,2)=watermarked_image_g;
imwrite(original_image,'withwatermarked_image.bmp','bmp');
%��������ʱ��
elapsed_time=cputime-start_time

%����psnr
PSNR=psnr(cover_object,withmark_image);

% ��ʾˮӡ��Ƕ��ˮӡͼ����ԭʼͼ��
figure(1)
subplot(2,2,2)
imshow(message);
title('ԭʼˮӡ');
subplot(2,2,3);
imshow(tempImg);
title('����ˮӡ');
subplot(2,2,4)
imshow(original_image);
name='Ƕ��ˮӡͼ��';
title(strcat(num2str(name),'   k=',num2str(k),'   PSNR=',num2str(PSNR)));
subplot(2,2,1)
imshow(original_image)
title('ԭʼͼ��');



