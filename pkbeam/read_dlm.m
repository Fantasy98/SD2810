function dQip= read_dlm(Z)
% Input: Z: null space
    
    load("dlmqip");
    niter = length(Qdlm.ktab);
    for k = 1:niter
        dQip.Qtab(:,:,k) = Z' * Qdlm.Qtab(:,:,k) * Z;
        dQip.ktab(k) = Qdlm.ktab(k);
    end
    dQip.rho = Qdlm.rho;
    dQip.bref = Qdlm.bref;
end