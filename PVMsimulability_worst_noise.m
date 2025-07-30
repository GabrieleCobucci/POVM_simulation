function v = PVMsimulability_worst_noise(d,N,target_povm)

%-------------------------------------------------------------------------%
%This function evaluates the SDP for PVM-simulability under worst-case
%noise.
%
%Inputs:
% - d: dimension of the POVM;
% - N: number of outcomes;
% - target_povm(:,:,a): POVM;
%
%Output:
% - v: critical visibility for PVM-simulation under worst-case noise.
%-------------------------------------------------------------------------%


%Identity
id=eye(d);

%Rank vector
plist=genpart(N,d);
Nr=size(plist,1);

%SDP
cvx_begin sdp quiet
cvx_solver mosek
variable F(d,d,N,Nr) hermitian semidefinite
variable noise(d,d,N) hermitian semidefinite
variable v

for a = 1 : N
    M(:,:,a) = v*target_povm(:,:,a) + noise(:,:,a);
end

for r = 1 : Nr
    s = 0;
    for a = 1 : N
        s = s + F(:,:,a,r);
    end
    q(r) = trace(s)/d;
end


for r = 1 : Nr
    s = 0; 
    for a = 1 : N
        s = s + F(:,:,a,r);
        trace(F(:,:,a,r)) == plist(r,a)*q(r);
    end
    s == q(r)*eye(d);
end

for a = 1 : N
    s = 0;
    for r = 1 : Nr
        s = s + F(:,:,a,r);
    end
    s == M(:,:,a);
end

sumnoise = 0;
for a = 1 : N
    sumnoise = sumnoise + noise(:,:,a);
end

sumnoise == (1-v)*id;


maximise(v)

cvx_end

end