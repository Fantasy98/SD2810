function [geo] = new_wing(y, xle, xte, ca, xea, xcm, my, Jy, EI, GK, nelem)
% [geo] = new_wing(y, xle, xte, ca, xea, xcm, my, Jy, EI, GK, nelem)
%
% Defines geometric and structural properties for a wing model
% based on beam theory. This is the most general version which
% assumes that all properties vary along the span.
%
% Input arguments:
% y     : Vector of spanwise positions at which the remaining
%         properties are given. These do not necessarily coincide
%         with the beam FEM nodes.
% xle   : x-coordinates of the leading edge at y
% xte   : x-coordinates of the trailing edge at y
% ca    : aileron chord at y (ca(i) = 0 if none present at y(i))
% xea   : x-coordinate of the elastic axis at y
% xcm   : x-coordinate of the center of mass at y
% my    : mass per length at y [kg/m]
% Jy    : torsional inertia per length at y [kgm?/m]
% nelem : number of finite elements required
%
% Output struct:
% geo.y    : spanwise position of the finite element nodes (vector)
% geo.x    : chordwise position of the midchord points (vector)
% geo.b    : semichord at FE nodes (vector)
% geo.bref : reference semichord (scalar)
% geo.ba   : aileron semichord at FE nodes, zero if no aileron (vector)
% geo.xea  : chordwise position of elastic axis (vector)
% geo.xcm  : chordwise position of center of mass (vector)
% geo.my   : mass per length at nodes [kg/m] (vector)
% geo.jy   : torsional inertia per length at nodes [kgm] (vector)
% geo.ei   : bending stiffness at nodes [Nm^2] (vector)
% geo.gk   : torsional stiffness at nodes [Nm^2] (vector)
% geo.rhof : fluid density, default 1.225 [kg/m^3] (scalar)
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % first, compute planform values for spanwise node
  % positions and chord dimensions
  yt = linspace(min(y), max(y), nelem+1);
  xl = interp1(y, xle, yt, 'spline');
  xt = interp1(y, xte, yt, 'spline');
  % vb = 0.5*(xt - xl);

  % fill struct with planform definition
  geo.y = yt;
  geo.x = 0.5*(xl + xt);
  geo.b = 0.5*(xt - xl);

  % print area (checking)
  area = trapz(yt, xt-xl);
  fprintf(1, 'Wing area: %.3f m^2\n', area);

  % semichord used as reference value
  geo.bref = 0.5*area / (max(yt) - min(yt));

  % compute aileron semichord vector (no interpolation)
  geo.ba = 0.5 * interp1(y, ca, yt, 'nearest');

  % interpolate mass properties
  geo.xea = interp1(y, xea, yt, 'spline');
  geo.xcm = interp1(y, xcm, yt, 'spline');
  geo.my = interp1(y, my, yt, 'spline');
  geo.jy = interp1(y, Jy, yt, 'spline');

  % the default wing does not have any lumped masses
  % Which could be changed according to the case?
  geo.cmi = [];
  geo.cdm = [];
  geo.cmx = [];
  geo.cmy = [];
  geo.cdj = {};

  % print total mass and inertia (checking)
  mass = trapz(yt, geo.my);
  tinert = trapz(yt, geo.jy);
  fprintf(1, 'Wing mass: %.3f kg\n', mass);
  fprintf(1, 'Torsional inertia: %.3f kgm^2\n', tinert);

  % interpolate structural properties
  geo.ei = interp1(y, EI, yt, 'spline');
  geo.gk = interp1(y, GK, yt, 'spline');

end

