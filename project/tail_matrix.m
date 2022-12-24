function [Qtail] = tail_matrix(k,nelem)
# Generate a 
    nnodes = nelem + 1;
    ndof = 3* nnodes;
    

    et1 = [1 0 0 1 0 0]';
    et3 = [0 0 1 0 0 1]';

    # Fixed parameter of HZT in report
    Stail= 0.86;
    le = 2.15;
    bt = 0.5 * Stail/le;
    # relative loc of elastic axis
    xa = 0.05;
    
    # The reduced frequency is 
    ke = k*bt/0.32; # 0.32 is semi-chord of wing (b)
    Qetail = beam_amatrix(ke,le,bt,xa);
    # Distance between AC of tail and AC of wing
    xt  = 2.64;

     # Tansformation matrix T2 
    T2 = [  1 -le/2 -xt
            0    1   0  
            0    0   1
            1  le/2  -xt 
            0    1   0 
            0    0   1];
    
    #Find location of fuselarge nodes
    Qtail = zeros(ndof,ndof);
    ivec = [3*fix(nnodes/2)+1 : 3*fix(nnodes/2)+3];

    # Give the contribution of tail 
    # Where T1 could be T2'
    Qtail(ivec,ivec) = T2' * Qetail * T2;


end