clear all;
 I0 = im2double(imread('frame11.png'));
 I1 = im2double(imread('frame10.png'));

%nomarlize the image.[0,1]
I0 = single((I0-min(I0(:))) ./ (max(I0(:)) - min(I0(:))));
I1 = single((I1-min(I1(:))) ./ (max(I1(:)) - min(I1(:))));
figure,
imshowpair(I0,I1);
title('before optical flow');

[M,N] = size(I0);
levels = [4,2,1];
um = zeros(2,2,'single'); vm = um;
pli = zeros(2,2,'single');

for gp = levels
    
    %smooth the image of current level.
    fs = fspecial('gaussian',[5 5],0.8);
    
    %resize the image to current level.
    wI0 = imresize( conv2(double(I0),fs,'same'),size(I0)/gp,'bilinear');
    wI1 = imresize( conv2(double(I1),fs,'same'),size(I0)/gp,'bilinear');
    [M,N] = size(wI0);
    %init parameter P.
    p11 = imresize(pli,size(I0)/gp,'bilinear');
    p12 = p11; p21 = p11; p22 = p11;
    
    [I1x,I1y] = gradient(wI1); %center gradient of vol2f
    
    [ um,vm ] = resizeFlow( um,vm,size(I0)./gp );
   
    for warped = 1:10 %current warping.
        D=zeros(M,N,2);
        D(:,:,1) = um; D(:,:,2) = vm;
        Iw1 = imwarp(wI1,D); % the effect of bound.
        
        Iw1x = imwarp(I1x,D);
        Iw1y = imwarp(I1y,D);
 %       [P,Iw1x,Iw1y] = sum_of_difference(Iw1,wI0);
        grad = Iw1x.*Iw1x+Iw1y.*Iw1y;
        rho_c = (Iw1 - Iw1x.*um-Iw1y.*vm-wI0);
        [um,vm,p11,p12,p21,p22] = tvl1_optimization(um,vm,grad,rho_c,Iw1x,Iw1y,p11,p12,p21,p22);    
    end
    
end
D = zeros(M,N,2);D(:,:,1)=um; D(:,:,2)=vm;
warped_I1 = imwarp(I1,D);
figure,
imshowpair(I0,warped_I1);
title('after optical flow');
figure,
imshow(flowToColor(D));
title('optical flow');