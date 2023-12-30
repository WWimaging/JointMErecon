function b=adj_forward_lr(x,s)
% d: k-space data
% s: sensitivity map

[Nx,Ny,~]=size(s);

a=ifft2(x);
b=a(1:Nx,1:Ny,:,:);




