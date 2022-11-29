function [Km,Mm,Zm,mQip]=ReduceDim(M,K,Qip,nmode)

% Input: 
%% M: Mass matrix with constrain
%% K: Stiffness matrix with constrain
%% Qip: Aerodynamic force with reduced frequency k 
%% nmode: Number of mode to implement 

%Output:
%% Km: reduced K matrix which shape=nmode x nmode
%% Mm: reduced M matrix which shape=nmode x nmode
%% Zm: projection 
%% mQip: reduced Qip 


% M: Mass matrix with constrain
    [Vre, D] = eig(K,M);
% Just use first 6 modes

    Zm = Vre(:,1:nmode);

% Retrieve reduced matrix
    Km = Zm' * K * Zm;
    Mm = Zm' * M * Zm;

% A class of aerodynamic force mQip
% fieldnames(mQip) to check the keys in class
    num_Q = size(Qip.Qtab,3);
    for i = 1:num_Q
        mQip.Qtab(:,:,i) = Zm' * Qip.Qtab(:,:,i) * Zm;
        mQip.ktab(i) = Qip.ktab(i);
    end
    mQip.bref = Qip.bref;
    mQip.rho = Qip.rho;

end