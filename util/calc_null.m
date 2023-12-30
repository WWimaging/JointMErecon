function [out]=calc_null(x,l_r,nc,nseg,pz)


[u,~,~] = svd(x,'econ');
un=u(:,(l_r+1):end);
un=reshape(un,pz,2*nc*nseg,[]);
un2=un(:,1:2:end,:)+1i*un(:,2:2:end,:);
un2=reshape(un2,pz*nc*nseg,[]);
out=un2.';

end