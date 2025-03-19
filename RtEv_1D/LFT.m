function [M11,M12,M21,M22]=LFT(n_psi,n_p,Nd)
    M11=zeros(n_p*(4*n_psi+Nd),n_p*(4*n_psi+Nd));
    M12=blkdiag(kron(ones(n_p,1),eye(n_psi)), ...
        kron(ones(n_p,1),eye(n_psi)), ...
        kron(ones(n_p,1),eye(n_psi)), ...
        kron(ones(n_p,1),eye(n_psi)), ...
        kron(ones(n_p,1),eye(Nd)));
    M21=blkdiag([zeros(n_psi,n_psi*n_p);eye(n_psi*n_p)], ...
        [zeros(n_psi,n_psi*n_p);eye(n_psi*n_p)], ...
        [zeros(n_psi,n_psi*n_p);eye(n_psi*n_p)], ...
        [zeros(n_psi,n_psi*n_p);eye(n_psi*n_p)], ...
        [zeros(Nd,Nd*n_p);eye(Nd*n_p)]);
    M22=blkdiag([eye(n_psi);zeros(n_p*n_psi,n_psi)], ...
        [eye(n_psi);zeros(n_p*n_psi,n_psi)], ...
        [eye(n_psi);zeros(n_p*n_psi,n_psi)], ...
        [eye(n_psi);zeros(n_p*n_psi,n_psi)], ...
        [eye(Nd);zeros(n_p*Nd,Nd)]);
end