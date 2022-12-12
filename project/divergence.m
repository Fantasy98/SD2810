function [udiv,zdiv] = divergence(K,Qip)
% [udiv,zdiv] = divergence(K,Qip)
%
% Compute the divergence speed and the corresponding
% deformation shape.
%
% K:     stiffness matrix
% Qip:   aerodynamic load interpolation struct
% udiv:  divergence speed
% zdiv:  divergence deformation shape
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % solve numerically well-posed divergence eigenvalue problem
  % for the critical dynamic pressure qdiv
  A = Qip.Qtab(:,:,1);
  [Modes,Pressure] = eig(A,K);
  [Lambda_max,indx] = max(abs(diag(Pressure)));
  qDiv = 1/Lambda_max;

  fprintf("\nThe divergence pressure is %.2f pa \n",qDiv);
  air_rho = Qip.rho;
  UDiv = sqrt(2*qDiv/air_rho);


  % determine critical speed from dynamic pressure
  udiv = UDiv;

  % corresponding eigenvector
  zdiv = Modes(:,indx);

end
