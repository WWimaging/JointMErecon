%% Example script: joint multi-echo reconstrution with structural low rank matrix completion
% The images of multiple echoes are jointly recovered
% wenchuan wu (wenchuan.wu@ndcn.ox.ac.uk)

addpath util
load data.mat
% Input
% rEPIall: the multi-echo k-space data [Nx Ny Nseg]
% M_sampling: sampling mask for all echoes [Nx Ny Nseg]
% senseMap: sensitivity maps [Nx Ny Nch]

l_rank=40;%rank parameter
l_lambda=1e-3;%regularisation parameter for SLR constraint
l_kernel=2;%kernel size for structured matrix formulation
n_iter=15;%number of iterations

[Nx,Ny,Nch,Nseg] = size(rEPIall);

% Direct ifft
im_ifftrecon=ifft2c(rEPIall);
im_ifftrecon_sos=(sum(abs(im_ifftrecon.^2),3)).^(1/2);
figure;montage(abs(im_ifftrecon_sos),"size",[1,3],"DisplayRange",[0 1e-3])
title("Direct ifft reconstruction of 3 rEPI echoes (from left to right, echo 1,2,3)");
set(gca,'FontName','Arial','FontSize',15)

% Sense recon
op_sense_me = @(x)adj_forward_sense_me(forward_sense_me(x,senseMap,M_sampling),senseMap,M_sampling);
b=adj_forward_sense_me(rEPIall,senseMap,M_sampling);
sense_recon = pcg(op_sense_me,b,1e-7,200);
sense_recon_reshape = reshape(sense_recon,[Nx,Ny,Nseg]);

% Joint recon
M_full=ones(size(M_sampling));
im_merecon = sense_recon(:);% use SENSE recon as the initial
pz=(2*l_kernel+1)^2;
for ii = 1 : n_iter
    ii
    im_merecon_curr=im_merecon;
    pxx=forward_sense_me(im_merecon,senseMap,M_full);
    M_SLR=mat_SLR(pxx,Nx,Ny,Nch*Nseg,l_kernel); 
    M_nu=calc_null(M_SLR,l_rank,Nch,Nseg,pz);
    [nfs,nfc]=ffilter(M_nu,Nx,Ny,Nch,Nseg,l_kernel);
    op_joint = @(x)forward_recon(x,senseMap,M_sampling,nfc,nfs,l_kernel,l_lambda);
    [im_merecon] = pcg(op_joint, b);
    thr = (norm(im_merecon_curr-im_merecon)/norm(im_merecon));
    if thr < 1e-6
        break;
    end
end

im_merecon_reshape=reshape(im_merecon,[Nx,Ny,Nseg]);

figure;montage(abs(sense_recon_reshape),"size",[1,3],"DisplayRange",[0 1e-3])
title("SENSE reconstruction of 3 rEPI echoes");
set(gca,'FontName','Arial','FontSize',15)
figure;montage(abs(im_merecon_reshape),"size",[1,3],"DisplayRange",[0 1e-3])
title("Joint reconstruction of 3 rEPI echoes");
set(gca,'FontName','Arial','FontSize',15)

