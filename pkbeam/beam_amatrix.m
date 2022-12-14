% Function for computing a single aerodynamics force matrix 
% Which will be implemented in Final Project
function [Ae] = beam_amatrix(k, le, b, xa)
    % [Ae] = beam_amatrix(k, le, b, xa)
    %
    % Construct unsteady aerodynamic load matrix Ae for a
    % beam element with length le and semichord b. 
    # xa is the elastic axis location w.r.t mid chord
    % local position of the elastic axis relative the the
    % midchord (xea - xmid).
    %
    % (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>
    
      % nondimensional distance xea - midchord
      a = xa/b;
    
      % Theodorsen's function
      if k == 0
        Ck = 1;
      else
        H0 = besselh(0,2,k);
        H1 = besselh(1,2,k);
        Ck = H1/(H1+i*H0);
      end
    
      % Theodorsen's frequency-domain aerodynamic coefficients
      Lh = -pi*((i*k)^2+2*Ck*(i*k));
      La = -pi*(a*(i*k)^2-(1+2*(1/2-a)*Ck)*(i*k)-2*Ck);
      Mh = -(pi/2)*(a*(i*k)^2+2*(1/2+a)*Ck*(i*k));
      Ma = -(pi/2)*((1/8+a^2)*(i*k)^2+((1/2-a)-2*(1/2+a)*(1/2-a)*Ck)*(i*k) -2*(1/2+a)*Ck);
    
      Ae = zeros(6,6);
      le2 = le^2;
      le3 = le^3;
    
      Ae(1,1) = (13/35)*(2*b)*(Lh/b)*le;
      Ae(1,2) = (11/210)*(2*b)*(Lh/b)*le2;
      Ae(1,3) = (7/20)*(2*b)*La*le;
      Ae(1,4) = (9/70)*(2*b)*(Lh/b)*le;
      Ae(1,5) = -(13/420)*(2*b)*(Lh/b)*le2;
      Ae(1,6) = (3/20)*(2*b)*La*le;
    
      Ae(2,1) = (11/210)*(2*b)*(Lh/b)*le2;
      Ae(2,2) = (1/105)*(2*b)*(Lh/b)*le3;
      Ae(2,3) = (1/20)*(2*b)*La*le2;
      Ae(2,4) = (13/420)*(2*b)*(Lh/b)*le2;
      Ae(2,5) = -(1/140)*(2*b)*(Lh/b)*le3;
      Ae(2,6) = (1/30)*(2*b)*La*le2;
    
      Ae(3,1) = (7/20)*(2*b)^2*(Mh/b)*le;
      Ae(3,2) = (1/20)*(2*b)^2*(Mh/b)*le2;
      Ae(3,3) = (1/3)*(2*b)^2*Ma*le;
      Ae(3,4) = (3/20)*(2*b)^2*(Mh/b)*le;
      Ae(3,5) = -(1/30)*(2*b)^2*(Mh/b)*le2;
      Ae(3,6) = (1/6)*(2*b)^2*Ma*le;
    
      Ae(4,1) = (9/70)*(2*b)*(Lh/b)*le;
      Ae(4,2) = (13/420)*(2*b)*(Lh/b)*le2;
      Ae(4,3) = (3/20)*(2*b)*La*le;
      Ae(4,4) = (13/35)*(2*b)*(Lh/b)*le;
      Ae(4,5) = -(11/210)*(2*b)*(Lh/b)*le2;
      Ae(4,6) = (7/20)*(2*b)*La*le;
    
      Ae(5,1) = -(13/420)*(2*b)*(Lh/b)*le2;
      Ae(5,2) = -(1/140)*(2*b)*(Lh/b)*le3;
      Ae(5,3) = -(1/30)*(2*b)*La*le2;
      Ae(5,4) = -(11/210)*(2*b)*(Lh/b)*le2;
      Ae(5,5) = (1/105)*(2*b)*(Lh/b)*le3;
      Ae(5,6) = -(1/20)*(2*b)*La*le2;
    
      Ae(6,1) = (3/20)*(2*b)^2*(Mh/b)*le;
      Ae(6,2) = (1/30)*(2*b)^2*(Mh/b)*le2;
      Ae(6,3) = (1/6)*(2*b)^2*Ma*le;
      Ae(6,4) = (7/20)*(2*b)^2*(Mh/b)*le;
      Ae(6,5) = -(1/20)*(2*b)^2*(Mh/b)*le2;
      Ae(6,6) = (1/3)*(2*b)^2*Ma*le;
end
    