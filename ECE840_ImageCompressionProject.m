%% Q5 (TV color image compression)*
clf
clc
clear all
 
% Transform the RGB image 'konark.jpg' given to you, into the YIQ space, which is used in NTSC television systems as:
 
Convert = [+0.299 +0.587 +0.114;
           +0.596 -0.274 -0.321;
           +0.211 -0.523 +0.311];
 
% [Y ; I; Q] = Matrix [R; G; B]
 
         
% Note that 0<R,G,B<1. (Wikipedia has more information: http://en.wikipedia.org/wiki/YIQ ).      
% Scale the YIQ components so that they lie within [0 1]. Then apply SVD based image compression.
 
% (i) Method-1 (all channels equally compressed): Use the first 24 singular values for Y, I, and Q.
 
subplot(4,2,1)
Cimg = imread('konark.jpg');
image(Cimg); axis off;
title('RGB Image')
Cimg = double(Cimg);
Cimg = Cimg/255;
 
YIQImg=zeros(size(Cimg));
for i=1:size(Cimg,1)
    for j=1:size(Cimg,2)
        YIQImg(i,j,1)=0.299*Cimg(i,j,1)+0.587*Cimg(i,j,2)+0.114*Cimg(i,j,3);
        YIQImg(i,j,2)=0.596*Cimg(i,j,1)-0.274*Cimg(i,j,2)-0.321*Cimg(i,j,3);
        YIQImg(i,j,3)=0.211*Cimg(i,j,1)-0.523*Cimg(i,j,2)+0.311*Cimg(i,j,3);
    end
end
 
YIQImg = double(YIQImg);
YIQImgDisp = YIQImg;
 
YImg = YIQImg(:,:,1);
IImg = YIQImg(:,:,2);
QImg = YIQImg(:,:,3);


% Y Image only
 
[U S V] = svd(YImg);
N=24;
ImgYpart1 = zeros(size(YImg));
ImgYpart1 = U(:,1:N)*S(1:N,1:N)*V(:,1:N)';
[M N] = size(YImg);
ImgYStart = zeros(M, N, 3);
ImgYStart(:,:,1) = YImg; ImgYStart(:,:,2) = YImg; ImgYStart(:,:,3) = YImg;
ImgYCompress = zeros(M, N, 3);
ImgYCompress(:,:,1) = ImgYpart1; ImgYCompress(:,:,2) = ImgYpart1; ImgYCompress(:,:,3) = ImgYpart1;
 
 
subplot(4,2,3); 
image(ImgYStart); axis off;
title('Y Image Start')
 
subplot(4,2,4); 
image(ImgYCompress); axis off;
title('Y Image Compressed Part 1')


% I Image Only
 
[U S V] = svd(IImg);
 
N=24;
ImgIpart1 = zeros(size(IImg));
ImgIpart1 = U(:,1:N)*S(1:N,1:N)*V(:,1:N)';
[M N] = size(IImg);
ImgIStart = zeros(M, N, 3);
ImgIStart(:,:,1) = IImg; ImgIStart(:,:,2) = IImg; ImgIStart(:,:,3) = IImg;
ImgICompress = zeros(M, N, 3);
ImgICompress(:,:,1) = ImgIpart1; ImgICompress(:,:,2) = ImgIpart1; ImgICompress(:,:,3) = ImgIpart1;
 
 
ImgIStart = abs(ImgIStart);
subplot(4,2,5); 
image(ImgIStart); axis off;
title('I Image Start')
 
ImgICompress = abs(ImgICompress);
subplot(4,2,6); 
image(ImgICompress); axis off;
title('I Image Compressed Part 1')

% Q Image Only
 
[U S V] = svd(QImg);
 
N=24;
ImgQpart1 = zeros(size(QImg));
ImgQpart1 = U(:,1:N)*S(1:N,1:N)*V(:,1:N)';
[M N] = size(IImg);
ImgQStart = zeros(M, N, 3);
ImgQStart(:,:,1) = QImg; ImgQStart(:,:,2) = QImg; ImgQStart(:,:,3) = QImg;
ImgQCompress = zeros(M, N, 3);
ImgQCompress(:,:,1) = ImgQpart1; ImgQCompress(:,:,2) = ImgQpart1; ImgQCompress(:,:,3) = ImgQpart1;
 
ImgQStart = abs(ImgQStart)
subplot(4,2,7); 
image(ImgQStart); axis off;
title('Q Image Start')
 
ImgQCompress = abs(ImgQCompress);
subplot(4,2,8); 
image(ImgQCompress); axis off;
title('Q Image Compressed Part 1')
 
% Compressed color
%convert back to RGB
 
r1 = 1.000*ImgYpart1+0.956*ImgIpart1+0.621*ImgQpart1;
g1 = 1.000*ImgYpart1-0.272*ImgIpart1-0.647*ImgQpart1;
b1 = 1.000*ImgYpart1-1.106*ImgIpart1+1.703*ImgQpart1;
RGBCompressed = zeros(M, N, 3); 
RGBCompressed(:,:,1) = r1; RGBCompressed(:,:,2) = g1; RGBCompressed(:,:,3) = b1;
RGBCompressed = abs(RGBCompressed);
 
subplot(4,2,2);
image(RGBCompressed); axis off;
title('RGB Image Compressed Part 1')

% (ii) Method-2 (chrominance compressed more): Use the first 52 singular values for Y and the first 10 singular values for I and Q each.
 
% Y Image only
figure%('Part 2 Compression')
 
subplot(4,2,1)
image(Cimg); axis off;
title('RGB Image')
 
 
[U S V] = svd(YImg);
N=52;
ImgYpart2 = zeros(size(YImg));
ImgYpart2 = U(:,1:N)*S(1:N,1:N)*V(:,1:N)';
[M N] = size(YImg);
ImgYStart = zeros(M, N, 3);
ImgYStart(:,:,1) = YImg; ImgYStart(:,:,2) = YImg; ImgYStart(:,:,3) = YImg;
ImgYCompress = zeros(M, N, 3);
ImgYCompress(:,:,1) = ImgYpart2; ImgYCompress(:,:,2) = ImgYpart2; ImgYCompress(:,:,3) = ImgYpart2;
 
subplot(4,2,3); 
image(ImgYStart); axis off;
title('Y Image Start')
 
subplot(4,2,4); 
image(ImgYCompress); axis off;
title('Y Image Compressed Part 2')
 
% I Image Only
 
[U S V] = svd(IImg);
 
N=10;
ImgIpart2 = zeros(size(IImg));
ImgIpart2 = U(:,1:N)*S(1:N,1:N)*V(:,1:N)';
[M N] = size(IImg);
ImgIStart = zeros(M, N, 3);
ImgIStart(:,:,1) = IImg; ImgIStart(:,:,2) = IImg; ImgIStart(:,:,3) = IImg;
ImgICompress = zeros(M, N, 3);
ImgICompress(:,:,1) = ImgIpart2; ImgICompress(:,:,2) = ImgIpart2; ImgICompress(:,:,3) = ImgIpart2;
 
ImgIStart = abs(ImgIStart);
subplot(4,2,5); 
image(ImgIStart); axis off;
title('I Image Start')
 
ImgICompress = abs(ImgICompress);
subplot(4,2,6); 
image(ImgICompress); axis off;
title('I Image Compressed Part 2')
 
% Q Image Only
 
[U S V] = svd(QImg);
 
N=10;
ImgQpart2 = zeros(size(QImg));
ImgQpart2 = U(:,1:N)*S(1:N,1:N)*V(:,1:N)';
[M N] = size(IImg);
ImgQStart = zeros(M, N, 3);
ImgQStart(:,:,1) = QImg; ImgQStart(:,:,2) = QImg; ImgQStart(:,:,3) = QImg;
ImgQCompress = zeros(M, N, 3);
ImgQCompress(:,:,1) = ImgQpart2; ImgQCompress(:,:,2) = ImgQpart2; ImgQCompress(:,:,3) = ImgQpart2;
 
ImgQStart = abs(ImgQStart);
subplot(4,2,7); 
image(ImgQStart); axis off;
title('Q Image Start')
 
ImgQCompress = abs(ImgQCompress);
subplot(4,2,8); 
image(ImgQCompress); axis off;
title('Q Image Compressed Part 2')

% RGB Compressed
 
subplot(4,2,2); 
RGBCompressed2 = zeros(M,N,3)
 
r2 = 1.000*ImgYpart2+0.956*ImgIpart2+0.621*ImgQpart2;
g2 = 1.000*ImgYpart2-0.272*ImgIpart2-0.647*ImgQpart2;
b2 = 1.000*ImgYpart2-1.106*ImgIpart2+1.703*ImgQpart2;
 
YIQCompressed2 = zeros(M, N, 3); 
YIQCompressed2(:,:,1) = r2; YIQCompressed2(:,:,2) = g2; YIQCompressed2(:,:,3) = b2;
YIQCompressed2 = abs(YIQCompressed2);
 
RGBCompressed2=abs(YIQCompressed2)
 
image(RGBCompressed2); axis off;
title('RGB Image Compressed Part 2')
 
% Figure of two compressed images
 
figure
subplot(2,1,1);
image(RGBCompressed); axis off;
title('RGB Image Compressed Part 1')
subplot(2,1,2);
image(RGBCompressed2); axis off;
title('RGB Image Compressed Part 2')


