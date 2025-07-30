function res = most_non_proj(d,N,Mstart)

%-------------------------------------------------------------------------%
%This function computes the numerical search for the most non-projective
%POVM under depolarising noise.
%
%Inputs:
% - d: dimension of the POVM;
% - N: number of outcomes;
% - Mstart(:,:,a): starting POVM for the numerical search;
%
%Output:
% - res: vector of critical visibilities at every step of the numerical
% search.
%-------------------------------------------------------------------------%

%Identity
id = eye(d);

v1 = 1;
v2 = 0;
res = [];

for k = 1 : N
    M(:,:,k) = Mstart(:,:,k);
end

while abs(v1-v2) > 1e-8

    v1 = v2;

    %First SDP: finding witness
    [v2,W] = PVMsimulability_dual(d,N,M);

    %Second SDP: maximum witness violation
    cvx_begin sdp quiet
    cvx_solver mosek
    variable M(d,d,N) hermitian semidefinite

    s = 0;
    obj = 1;
    for a = 1 : N
        s = s + M(:,:,a);
        obj = obj + trace(M(:,:,a)*W(:,:,a));
    end

    s == id;

    minimise(obj)
    cvx_end

    res = [res, v2]
end



end