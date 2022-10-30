
% lwmain.m
%
% main program for aeroelastic analysis
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

clear all;

% setup geometry and structural properties
% number of finite elements requested should be a multiple of 3
nelem = 24;
nnodes = nelem + 1;

% lab wing dimensions and properties
l =
b =
ba =
mhinge =
t =
rhop =
E =
G =

% definition matrix for discrete point masses to attach
npmass = ;
dpm = zeros(npmass,3);
dpm(1,:) =
dpm(2,:) =
% ....

% set up linear constraints for clamped wing root
B =

% retrieve system matrices
[M,K,Z,Qip,f,CRv,CRd] = labwing(B, l, b, t, ba, mhinge, rhop, E, G, nelem, dpm);

% print free vibration frequencies


% show modeshapes
plotmode( ??? );

% stop here until the rest is implemented
return;

% compute divergence speed
[udiv,zdiv] = divergence(K, Qip);
fprintf(1,'Divergence speed: %.2f m/s \n', udiv);

% compute reversal speed
[urev,zrev] = reversal(K, Qip, f, CRv, CRd);
fprintf(1,'Reversal speed: %.2f m/s \n', urev);

% compute flutter speed
[ucrit, pcrit, zcrit] = flutter(M,K,Qip);
fcrit =
fprintf(1,'Flutter speed: %.2f m/s \n',ucrit);
fprintf(1,'Frequency of the critical mode: %.2f Hz \n',fcrit);

% look at flutter mode shape
vismode( ??? );
