function [Q0,vtot,vhat] = flight_load(K,M,B,Qip,nelem,u,nz,k)
# function to compute the deformation of wing, including rigid body: alpha & delta and elastic deformatipn
# Inputs:
%  K: unconstrained stiffness matrix
%  M: unconstrained mass matrix
%  B: Constrain matrix
%  Qtab: aerodynamic force of wing based on strip theory 
%  nelem: number of elements
%  nz : load factor 
%  u: flight speed
%  k: reduced frequency





    nnodes = nelem + 1;
    ndof = 3* nnodes;

# Steady aerodynamic force
    Q0 = Qip.Qtab(:,:,1);


# indices vetors
    e1 = zeros(ndof,1);
    e1(1:3:end) = 1;
    e3 = zeros(ndof,1);
    e3(3:3:end) = 1;
    
    Z = null(B);
    # Retreive total mass
    mtot = e1' * M * e1;
    # Dynamic pressure
    q = 0.5*Qip.rho*u*u;  g = 9.81;
    

    #### Install Horizontal tail aerodynamic matrix
    # Assume only ONE element == 2 nodes
    et1 = [1 0 0 1 0 0];
    et3 = [0 0 1 0 0 1];

    # Fixed parameter of HZT in report
    Stail= 0.86;
    le = 2.15;
    bt = 0.5 * Stail/le;
    # relative loc of elastic axis
    xa = 0.05;

    Qetail = beam_amatrix(k,le,bt,xa);
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
    Qtail = zeros(size(Q0));
    ivec = [3*fix(nnodes/2)+1 : 3*fix(nnodes/2)+3];

    # Give the contribution of tail 
    # Where T1 could be T2'
    Qtail(ivec,ivec) = T2' * Qetail * T2;
    # Assemble aerodynamic matrix with tail   
    Q0 = Q0 + Qtail;

    
    ######### Solve the inertia relief 

    LHS = [ Z' *(K-q*Q0)*Z  -q*Z'*Q0*e3
            q*e1'*Q0*Z       q*e1'*Q0*e3];

    % Right hand side

    RHS = [-nz*g*Z'*M*e1 
            g*e1'*M*e1];

    x = LHS\RHS;

    # Extract elastic deformation 
    vhat = x(1:end-1);
    # Extract angle of attack as rigid body motion
    alpha0 = x(end);
    # Assemble the deformation vector 
    vtot = Z*vhat + alpha0 * e3;

    # Check residual of solution by the equlibrium
    res = K * vtot + M * nz * g* e1 - q * Q0 * vtot;
    fprintf("Computation residual =  %E \n",min(abs(res)))







end