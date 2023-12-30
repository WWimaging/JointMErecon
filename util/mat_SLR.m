function b = mat_SLR(a,Nx,Ny,Ns,l_k)
% Structural low rank matrix formulation - the 'S' matrix
% using 'S' type matrix in loraks
% a: input k space data
% Ns: Nc * Nseg
% l_k:kernel size

% wenchuan wu (wenchuan.wu@ndcn.ox.ac.uk)

a1=reshape(a,Nx,Ny,Ns);

pz=(l_k*2+1)^2;
nk=(Nx-l_k*2-1)*(Ny-l_k*2-1);%number of patches
m_slr=zeros(Ns,pz,nk);


cc =0;
for ii = l_k+2:(Nx-l_k)
    for jj=l_k+2:(Nx-l_k)
        cc=cc+1;
        mm=[];
        for i1=-l_k:l_k
            for j1=-l_k:l_k
                    mm=cat(2,mm,squeeze(a1(ii+i1,jj+j1,:)));
            end
        end
        m_slr(:,:,cc)=mm;

    end
end

m_slr_symm=flip(m_slr,3);

m_slr_S_tmp1=cat(4,(real(m_slr)-real(m_slr_symm)),(-imag(m_slr)+imag(m_slr_symm)));
m_slr_S_tmp2=cat(4,(imag(m_slr)+imag(m_slr_symm)),(real(m_slr)+real(m_slr_symm)));
b=cat(3,m_slr_S_tmp1,m_slr_S_tmp2);
b=reshape(permute(b,[2 4 1 3]),[Ns*pz*2,nk*2]);

end



