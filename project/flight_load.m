function [vtot] = flight_load(K,M,B,Qtab,nelem,nz)
# function to compute the deformation of wing, including rigid body: alpha & delta and elastic deformatipn
# Inputs:
%  K: unconstrained stiffness matrix
%  M: unconstrained mass matrix
%  B: Constrain matrix
%  Qtab: aerodynamic force of wing based on strip theory 
%  nelem: number of elements
%  nz : load factor 
%  u: flight speed

    nnode = nelem + 1;
    ndof = 3* nnode;

# Steady aerodynamic force
    Q0 = Qip.Qtab(:,:,1);
# indices vetors
    e1 = zeros(ndof,1);
    e1(1:3:end) = 1;
    e3 = zeros(ndof,1);
    e3(3:3:end) = 1;
    # Retreive total mass
    mtot = e1' * M * e1;
    q = 0.5*Qip.rho*u*u;  g = 9.81;
    # Compute the Lift when manuver since L = mg
    L = mtot*g;


















end