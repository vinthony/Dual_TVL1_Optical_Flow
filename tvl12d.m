clear all;clc;
 I0 = im2double(rgb2gray(imread('car1.jpg')));
 I1 = im2double(rgb2gray(imread('car2.jpg')));

%nomarlize the image.[0,1]
I0 = (I0-min(I0(:))) ./ (max(I0(:)) - min(I0(:)));
I1 = (I1-min(I1(:))) ./ (max(I1(:)) - min(I1(:)));
figure,
imshowpair(I0,I1);
title('before optical flow');

[M,N] = size(I0);
levels = [4,2,1];
um = zeros(2,2,'single'); vm = um;
pli = zeros(2,2,'single');
tau = 0.25;
theta = 0.1;
lambda =50;
taut = tau / theta;
OUT_ITERATION = 10;
INNER_ITERATION = 30;
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
        %figure,imshow(flowToColor(D));
        Iw1 = imwarp(wI1,D); % the effect of bound.
        Iw1x = imwarp(I1x,D);
        Iw1y = imwarp(I1y,D);
        grad = Iw1x.*Iw1x+Iw1y.*Iw1y;
        rho_c = (Iw1 - Iw1x.*um-Iw1y.*vm-wI0);
        
        n_inter=0;n_outer=0;error = 1000000;
        
        while error > epsilon*epsilon && n_outer < OUT_ITERATION
             n_outer = n_outer + 1;
%               um = medfilt2(um,[5,5]);
%               vm = medfilt2(vm,[5,5]);
            while error > epsilon*epsilon && n_inter < INNER_ITERATION
                n_inter = n_inter + 1;
                %TH
                [v1,v2] = ev2d(grad,rho_c,Iw1x,Iw1y,um,vm,lambda,theta); 
                div_p1 = mydiv2d(p11,p12);
                div_p2 = mydiv2d(p21,p22);
                um0 = um; vm0 = vm;
                um = v1 + theta .* div_p1;
                vm = v2 + theta .* div_p2;
                e = (um0-um).*(um0-um)+(vm0-vm).*(vm0-vm);
                error = sum(e(:))/12;
                errors(OUT_ITERATION*(n_outer-1)+n_inter,warped) = error;
                [umx,umy] = forward_gradient_2d(um);
                [vmx,vmy] = forward_gradient_2d(vm);
                g1 = sqrt(umx.*umx+umy.*umy);
                g2 = sqrt(vmx.*vmx+vmy.*vmy);
                ng1 = 1+taut .*g1;
                ng2 = 1+taut .*g2;
                p11 = (p11 + taut .* umx) ./ng1;
                p12 = (p12 + taut .* umy) ./ng1;
                p21 = (p21 + taut .* vmx) ./ng2;
                p22 = (p22 + taut .* vmy) ./ng2;
            end
        end
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