function plotmode_multi(v)
% plotmode(v)
%
% Plots nodal deflection and twist versus span.
%
% v :     deflection to plot, assuming beam DOFs [w, phi, theta]
%
% (c) 2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>

  % constants
  l = 1.6;

  % dimensions
  ndof = length(v);
  nnode = ndof/3;
  yp = linspace(0.0, l, nnode);

  % check dimension of v
  if numel(v) ~= ndof
    error('plotmode: Can only plot one modeshape at a time.');
  elseif mod(ndof,3) ~= 0
    error('plotmode: Number of DOFs is not a multiple of 3.');
  end

  % Extract bending and torsion dofs section by section
  for k = 1:nnode
    w(k) = v(1+3*(k-1)); % nodal deflection
    t(k) = 180/pi * v(3+3*(k-1)); % nodal twist
  end

  % distinct plot for no deformation
  if max(abs(w)) < 1e-9
    w = zeros(1,nnode);
  elseif max(abs(t)) < 1e-9
    t = zeros(1,nnode);
  end

  % Plot modeshape
  
  % subplot(4,1,[1 2])
  hold on
    plot(yp, t,"linewidth",7)
    % title('Deformation of elastic axis',"fontsize",25)
    ylab=ylabel('Deflection [m]',"fontsize",20);
    xlabe =xlabel('Span coordinate [m]',"fontsize",20);
    
    box
  hold off
  % subplot(4,1,[3 4])
  % hold on
  %   plot(yp, t,"linewidth",1.5)
    
  %   ylabel('Twist [deg]',"fontsize",20)
  %   box
  % hold off
  %for save figure
  %print -djpg Modeshape.jpg
