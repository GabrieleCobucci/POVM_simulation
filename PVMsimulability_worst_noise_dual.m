function [s,Gamma] = PVMsimulability_worst_noise_dual(d,N,target_povm)

%-------------------------------------------------------------------------%
%This function evaluates the dual SDP for PVM-simulability under worst-case
%noise
%
%Inputs:
% - d: dimension of the POVM;
% - N: number of outcomes;
% - target_povm(:,:,a): POVM;
%
%Output:
% - v: critical visibility for PVM-simulation under worst-case noise;
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
    variable B(d,d) hermitian
    
    for a = 1 : N
        - Gamma(:,:,a) -  B >= 0;
    end
    
    for r = 1 : Nr
        sumavar = 0;
        for a = 1 : N
            sumavar = sumavar - y(a,r)*plist(r,a);
        end
        sumavar - trace(Y(:,:,r)) == 0;
    end
    
    for r = 1 : Nr
        for a = 1 : N
            Gamma(:,:,a) - y(a,r)*id - Y(:,:,r) >= 0;
        end
    end
    
    s1 = 1;
    for a = 1 : N
        s1 = s1 + trace(Gamma(:,:,a)*target_povm(:,:,a));
    end
    
    s1 + trace(B) == 0;
    
    obj = -trace(B);
    
    minimise(real(obj))
    
    cvx_end
    
    s = cvx_optval;

end
