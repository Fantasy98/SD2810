function [Q,vtot,ve,alfa0,delta0] = flight_load(K,M,B,Qip,f,nelem,u,nz,k)
#Quasi-Steady aerodynamic 
#function to compute the deformation of wing, including rigid body: alpha & delta and elastic deformatipn
# Inputs:
%  K: unconstrained stiffness matrix
%  M: unconstrained mass matrix
%  B: Constrain matrix
%  Qtab: aerodynamic force of wing based on strip theory 
%  f: aileron load vector
%  nelem: number of elements
%  nz : load factor : for our case could be 5.3, 7, or 9 
%  u: flight speed
%  k: reduced frequency

# Outputs:
%   Q: Assembled aerodynamic matrix consider the horizontal tail
%   vtot: The total deflection vector size = [1, 3*(nelem+1)]
%   ve: elastic deflection vector size = [1, 3*(nelem+1)]
%   alfa0: reference angle of attack
%   delta0: aileron deflection

    nnodes = nelem + 1;
    ndof = 3* nnodes;
    
# Aerodynamic force based on reduced freq
    Q0 = ipolQk(Qip,k);
    

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
    
    Qtail = tail_matrix(k,nelem);
   
    Q = Q0 + Qtail;
    ######### Solve the inertia relief 

    LHS = [
            Z'*(K-q*Q)*Z    -q*Z'*Q*e3      -q*Z'*f
            
            q*e1'*Q*Z       q*e1'*Q*e3      q*e1'*f
            
            q*e3'*Q0*Z        q*e3'*Q0*e3      q*e3'*f
            % q*e3'*Q*Z        q*e3'*Q*e3      q*e3'*f
           
            ];

    % Right hand side

    RHS = [
            -Z'*M*e1*nz*g

            e1'*M*e1*g*nz
            
            0
            ];

    x = LHS\RHS;

    # Extract elastic deformation 
    vhat = x(1:end-2);
    # Extract angle of attack as rigid body motion
    alfa0 = x(end-1);
    delta0 = x(end);
    # Assemble the deformation vector 
    vtot = Z*vhat + alfa0*e3 + delta0*e3;
    ve = Z*x(1:end-2);
    
    # Check residual of solution by the equlibrium
    res = K * vtot + M * nz * g* e1 - q * Q0 * vtot;
    fprintf("Computation residual =  %E \n",min(abs(res)))







end