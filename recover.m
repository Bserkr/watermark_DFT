%��ȡ
clc
clear all;

% ���濪ʼʱ��
start_time=cputime;
iTimes=22;              %���Ҵ���
blocksize=8;            % ���ÿ�Ĵ�С
filter_m=[  1,1,1,1,1,1,1,1;    % �˲�����
            1,1,1,1,1,1,1,1;
            1,1,0,0,0,0,1,1;
            1,1,0,0,0,0,1,1;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;];
% ����Ƕ��ˮӡͼ��
withmark_1=imread('withwatermarked_image.bmp');
withmark=withmark_1(:,:,2);
imshow(withmark_1)
watermarked_image=double(withmark)/255;

% Ƕ��ˮӡͼ����������������
Mw=size(watermarked_image,1);	
Nw=size(watermarked_image,2);	        

% ����Ƕ����Ϣ��
max_message=Mw*Nw/(blocksize^2);

% ����ԭʼˮӡ
mark=imread('mark40.bmp');
mark=rgb2gray(mark);
orig_watermark=double(mark);

% ԭʼˮӡ���������������
Mo=size(orig_watermark,1);	
No=size(orig_watermark,2);	

%���������������״̬Ϊ1100
key=1100;
rand('state',key);

% ����α�������
pn_sequence_zero=round(2*(rand(1,sum(sum(filter_m)))-0.5));
pn_sequence_one=round(2*(rand(1,sum(sum(filter_m)))-0.5));              
% ��ͼ��ֿ�
x=1;
y=1;
h=waitbar(0,'��ȡˮӡ����ȴ�');
for (kk = 1:max_message)

    % ����Ҷ�任
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
   
    % ����sequence��pn_sequence_zero��pn_sequence_one�����ϵ��
   correlation_zero(kk)=corr2(pn_sequence_zero,sequence);
   correlation_one(kk)=corr2(pn_sequence_one,sequence);               
    
    % �ƶ�����һ��
    if (x+blocksize) >= Nw
        x=1;
        y=y+blocksize;
    else
        x=x+blocksize;
    end
    waitbar(kk/max_message,h);
end
close(h);

% ���correlation_zero>correlation_one����ômessage_vector=0����֮Ϊ1.
for (kk=1:Mo*No)
    if correlation_zero(kk)>correlation_one(kk)
        message_vector(kk)=0;
    else
        message_vector(kk)=1;
    end
end

% Arnold����
tempImg=reshape(message_vector(1:Mo*No),Mo,No);
message_arnold=tempImg;
for n=1:iTimes % ����
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

% ��������ʱ��
elapsed_time=cputime-start_time,
%����NC����һ�����ϵ����
NC=nc(message,orig_watermark)

% ��ʾ��ȡˮӡ��ԭʼˮӡ
figure(3)
subplot(1,2,1);
imshow(message,[]);
name='��ȡˮӡ';
title(strcat(num2str(name),'NC=',num2str(NC)));
subplot(1,2,2)
imshow(orig_watermark,[])
title('ԭʼˮӡ');