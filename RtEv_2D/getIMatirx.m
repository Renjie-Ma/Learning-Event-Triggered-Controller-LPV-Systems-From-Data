function [I1,I2,I3]=getIMatirx(n_psi,p_ch,Nd)
I1=[eye(n_psi);kron(p_ch,eye(n_psi))];
I2=[eye(Nd);kron(p_ch,eye(Nd))];
I3=[eye(n_psi);kron(p_ch,eye(n_psi));kron(p_ch,kron(p_ch,eye(n_psi)))];
end