function v = PVMsimulability(d,N,M)

%-------------------------------------------------------------------------%
%This function evaluates the SDP for PVM-simulability under depolarising
%noise.
%
%Inputs:
% - d: dimension of the POVM;
% - N: number of outcomes;
% - M(:,:,a): POVM;
%
%Output:
% - v: critical visibility for PVM-simulation under depolarising noise.
%-------------------------------------------------------------------------%

%Identity
id = eye(d);

%Rank vectors
plist = genpart(N,d);
Nr = size(plist,1);

%SDP
cvx_begin sdp quiet
cvx_solver mosek
variable F(d,d,N,Nr) hermitian semidefinite
variable v

v >= 0;

for r = 1 : Nr
    sum1 = 0;
    for a = 1 : N
        sum1 = sum1 + F(:,:,a,r);
    end
    q{r} = trace(sum1)/d;
    sum1 == q{r}*id;
end

for a = 1 : N
    for r = 1 : Nr
        trace(F(:,:,a,r)) == q{r}*plist(r,a);
    end
end

for a = 1 : N
    sumf = 0;
    for r = 1 : Nr
        sumf = sumf + F(:,:,a,r);
    end
    sumf == v*M(:,:,a) + (1-v)/d*trace(M(:,:,a))*id;
end

maximize(v)
cvx_end

end