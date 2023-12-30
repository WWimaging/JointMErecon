function y=forward_sense_me(x,s,m)
% x: image
% s: sensitivity map
% m: k-space sampling

[nx,ny,nc]=size(s);
nseg = length(x(:))/(nx*ny);
y=[];
ll=nx*ny;
for iseg = 1 : nseg

    x_t = x((ll*(iseg-1)+1):ll*iseg);
    x_t = reshape(x_t,[nx,ny]);
    x_t_mc = repmat(x_t,[1,1,nc]).*s;
    k_t_mc = fft2c(x_t_mc).*m(:,:,iseg);
    y=cat(1,y,k_t_mc(:));

end
end