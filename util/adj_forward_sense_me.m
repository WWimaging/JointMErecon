function y=adj_forward_sense_me(d,s,m)
% d: k-space data
% s: sensitivity map
% m: k-space sampling

[nx,ny,nc]=size(s);
nseg = length(d(:))/(nx*ny*nc);
y=[];
ll=nx*ny*nc;
for iseg = 1 : nseg

    d_t = d((ll*(iseg-1)+1):ll*iseg);
    d_t = reshape(d_t,[nx,ny,nc]).*m(:,:,iseg);
    d_t_sc = sum(ifft2c(d_t).*conj(s),3);
    y=cat(1,y,d_t_sc(:));

end

end