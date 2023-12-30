function y=forward_recon(x,s,m,n1,n2,l_k,l_l)

m_f = ones(size(m));

d_sense= adj_forward_sense_me(forward_sense_me(x,s,m),s,m);
d_lr= 2*adj_forward_sense_me(adj_forward_lr(forward_lr(x,s,m_f,n1,l_k),s)-adj_forward_lr(forward_lr(x,s,m_f,n2,l_k,1),s),s,m_f);

y=d_sense+l_l*d_lr;

end