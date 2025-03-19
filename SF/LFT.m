function [M11,M12,M21,M22]=LFT(n_x,n_p,Nd)
    M11=zeros(n_p*(4*n_x+Nd),n_p*(4*n_x+Nd));
    M12=blkdiag(kron(ones(n_p,1),eye(n_x)), ...
        kron(ones(n_p,1),eye(n_x)), ...
        kron(ones(n_p,1),eye(n_x)), ...
        kron(ones(n_p,1),eye(n_x)), ...
        kron(ones(n_p,1),eye(Nd)));
    M21=blkdiag([zeros(n_x,n_x*n_p);eye(n_x*n_p)], ...
        [zeros(n_x,n_x*n_p);eye(n_x*n_p)], ...
        [zeros(n_x,n_x*n_p);eye(n_x*n_p)], ...
        [zeros(n_x,n_x*n_p);eye(n_x*n_p)], ...
        [zeros(Nd,Nd*n_p);eye(Nd*n_p)]);
    M22=blkdiag([eye(n_x);zeros(n_p*n_x,n_x)], ...
        [eye(n_x);zeros(n_p*n_x,n_x)], ...
        [eye(n_x);zeros(n_p*n_x,n_x)], ...
        [eye(n_x);zeros(n_p*n_x,n_x)], ...
        [eye(Nd);zeros(n_p*Nd,Nd)]);
end