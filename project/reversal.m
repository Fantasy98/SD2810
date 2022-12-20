function [urev,zrev] = reversal(K,Q,f,CRv,CRd,Z)
% [urev,zrev] = reversal(dsy, geo)
%
% Compute the airpeed at aileron reversal and the deformation
% for unit flap deflection at reversal speed.
%
% geo:    wing geometry description
% dsy:    aeroleastic system
% urev:   critical airspeed at aileron reversal
% zrev:   wing deformation vector for unit delta at qrev
%
% (c) 2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

% Compute the aerodynamic force considering the control surface deflection
% Qip.Qtab(:,:,1) corresponding to the static state 
A = Q - f*(CRv./CRd);
% Solve for eigenvalue problem
[V,D] = eig(Z'*A*Z,Z'*K*Z);
% Find the maximum eigenvalue whose reverse is ths dynamic pressure
% and location as well
[lambda_rev,max_loc] = max(abs(diag(D)));
qrev = 1/lambda_rev;
% q = 0.5*rho*u^2
urev = sqrt(2*qrev/1.225);
zrev = V(:,max_loc);

end
