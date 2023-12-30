function [a,b]=ffilter(x,Nx,Ny,Nch,Nseg,l_k)

[n_f,~]=size(x);
xx=reshape(x,n_f,(l_k*2+1),(l_k*2+1),Nch*Nseg);
xx=permute(xx,[3,2,4,1]);

fil1=fft2(conj(xx),4*l_k+1,4*l_k+1);
cfil1=fil1;
ptch=permute(cfil1,[1,2,3,5,4]).*permute(fil1,[1,2,5,3,4]);
ptch=padarray(ifft2(sum(ptch,5)),[Nx-1-2*l_k,Ny-1-2*l_k],'post');
a=fft2(circshift(ptch,[-4*l_k -4*l_k]));


fil2=flip(flip(xx,1),2);
fil2=fft2(fil2,4*l_k+1,4*l_k+1);
ptch=permute(cfil1,[1,2,5,3,4]).*permute(fil2,[1,2,3,5,4]);
ptch=padarray(ifft2(sum(ptch,5)),[Nx-1-2*l_k,Ny-1-2*l_k],'post');
b=fft2(circshift(ptch,[-2*l_k -2*l_k]));


end


