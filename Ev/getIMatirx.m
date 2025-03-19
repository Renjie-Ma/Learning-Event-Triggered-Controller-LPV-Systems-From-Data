function [I1,I2,I3]=getIMatirx(n_x,p_ch,Nd)
I1=[eye(n_x);kron(p_ch,eye(n_x))];
I2=[eye(Nd);kron(p_ch,eye(Nd))];
I3=[eye(n_x);kron(p_ch,eye(n_x));kron(p_ch,kron(p_ch,eye(n_x)))];
end