function [umo,vmo,p11,p12,p21,p22] = tvl1_optimization(um,vm,grad,rho_c,Iw1x,Iw1y,p11,p12,p21,p22)
    epsilon = 1e-2;
    OUT_ITERATION = 10;
    INNER_ITERATION = 30;
    tau = 0.25;
    [M,N] = size(vm);
    sz = M*N;
    theta = 0.1;
    lambda =50;
    taut = tau / theta;
    n_inter=0;n_outer=0;error = 1000000;
    while error > epsilon*epsilon && n_outer < OUT_ITERATION
        n_outer = n_outer + 1;
        um = medfilt2(um,[5,5]);
        vm = medfilt2(vm,[5,5]);
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
            error = sum(e(:))/sz;
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
    umo = um; vmo = vm;
end