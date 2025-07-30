function [v,Gamma] = PVMsimulability_dual(d,N,M)

%-------------------------------------------------------------------------%
%This function evaluates the dual SDP for PVM-simulability under
%depolarising noise.
%
%Inputs:
% - d: dimension of the POVM;
% - N: number of outcomes;
% - M(:,:,a): POVM;
%
%Output:
% - v: critical visibility for PVM-simulation under depolarising noise;
% - Gamma: operators to construct the witness.
%-------------------------------------------------------------------------%


%Identity
id = eye(d);

%Rank vectors
plist = genpart(N,d);
Nr = size(plist,1);

%SDP
cvx_begin sdp quiet
cvx_solver mosek
variable Gamma(d,d,N) hermitian
variable Y(d,d,Nr) hermitian
variable y(N,Nr)

for r = 1 : Nr
    sumy = 0;
    for a = 1 : N
        Gamma(:,:,a) - Y(:,:,r) - y(a,r)*id >= 0;
        sumy = sumy + y(a,r)*plist(r,a);
    end
    - sumy - trace(Y(:,:,r)) == 0;
end

sum1 = 1;
sum2 = 0;
for a = 1 : N
    sum1 = sum1 + trace(Gamma(:,:,a)*M(:,:,a));
    sum2 = sum2 + 1/d*trace(M(:,:,a))*trace(Gamma(:,:,a));
end

sum1 == sum2;

minimise(real(sum1))
cvx_end

v = cvx_optval;

end