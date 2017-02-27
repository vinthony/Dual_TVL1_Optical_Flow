function [ um,vm ] = resizeFlow( um,vm,new_size )
%RESIZEFLOW Summary of this function goes here
%   Detailed explanation goes here
    % resize flow
    scaling = single([new_size(1)/size(um,1),new_size(2)/size(um,2)]);
    a = linspace(1,size(um,2),new_size(2));
    b = linspace(1,size(um,1),new_size(1));
    [xi,yi] = meshgrid(a,b);
    xi = single(xi);yi=single(yi);
    um = interp2(single(um).*scaling(2),xi,yi); um(isnan(um)) =0;
    vm = interp2(single(vm).*scaling(1),xi,yi); vm(isnan(vm)) =0;

end

