clear all;clc;
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
    % resize to current level. smooth
    fs = fspecial('gaussian',[5 5],0.8);
    wI0 = conv2(double(I0),fs,'same');
    wI1 = conv2(double(I1),fs,'same');
   
    wI0 = imresize(wI0,size(I0)/gp,'bilinear');
    wI1 = imresize(wI1,size(I0)/gp,'bilinear');
    p11 = imresize(pli,size(I0)/gp,'bilinear');
    p12 = p11; p21 = p11; p22 = p11;
    
    % resize flow
    new_size = size(I0)/gp;
    scaling = single([new_size(1)/size(um,1),new_size(2)/size(um,2)]);
    a = linspace(1,size(um,2),new_size(2));
    b = linspace(1,size(um,1),new_size(1));
    [xi,yi] = meshgrid(a,b);
    xi = single(xi);yi=single(yi);
    um = interp2(single(um).*scaling(2),xi,yi); um(isnan(um)) =0;
    vm = interp2(single(vm).*scaling(1),xi,yi); vm(isnan(vm)) =0;
    
    [I1x,I1y] = gradient(wI1); %center gradient of vol2f
    
    epsilon = 1e-2;
    MAX_ITERATION = 256;
    errors = zeros(MAX_ITERATION,10);
    for warped = 1:10 %current warping.
        D=zeros(new_size(1),new_size(2),2);
        D(:,:,1) = um; D(:,:,2) = vm;
        Iw1 = imwarp(wI1,D); % the effect of bound.
        Iw1x = imwarp(I1x,D);
        Iw1y = imwarp(I1y,D);
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