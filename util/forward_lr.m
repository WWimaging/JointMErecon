function b=forward_lr(x,s,m,nn,l_k,c_tag)
% x: k-space data
% s: sensitivity map
% m: k-space sampling
% nn:filter,C/S formulation
% l_k:kernel size
% c_tag:tag for C matrix

if nargin < 6
    c_tag=0;
end

[Nx,Ny,Nch]=size(s);
Nseg = size(m,3);


a=forward_sense_me(x,s,m);
b=padarray(reshape(a,[Nx Ny Nch*Nseg]),[2*l_k, 2*l_k], 'post');
if c_tag
b=squeeze(sum(nn.*repmat(conj(fft2(b)),[1 1 1 Nch*Nseg]),3));
else
b=squeeze(sum(nn.*repmat(fft2(b),[1 1 1 Nch*Nseg]),3));
end





