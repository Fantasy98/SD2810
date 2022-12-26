function [SF_n,SF_s] = safety_factor(v)

    l =12.1; % m

    c = 0.64;
# According to the technic report, the main wing is assumed as ellipse with chord = 500mm
    b = c/2 ; % m
# The alieron span is defined as c*E, E = 0.225 in W12C-JARA-0001.pdf P
    t = 6e-3;
    G = 8600E6;
    E = 23.9E9;
    
    Ixx = 1.6117e+05 * 10^(-12);
    
    ndof = length(v);
    nnode = fix(ndof/3);
    nshape = [1, nnode];
    nelem = nnode-1;
    # We need elemental deformation! See the equation (59) in the LectureNote


   
    % Extract bending and torsion dofs section by section
    w = zeros(nshape);
    theta = zeros(nshape);
    for k = 1:nnode
        w(k) = v(1+3*(k-1)); % nodal deflection
        w1(k) = v(2+3*(k-1)); % nodal deflection curvature

        
        theta(k) = 180/pi * v(3+3*(k-1)); % nodal twist
    end
    yp = linspace(0.0, l, nnode);
    ye = linspace(0,l,nelem);

    # Interpolation of w'' 
    le = yp(2)-yp(1);
    w2 = diff(w1)/le;
    

    t = c * 0.01;
    % t = 0.11;
    Mx = E * Ixx .* w2';
    % Z = Mx /(0.5*t) ;
    sigma = 0.5*t*Mx /Ixx;

    K = 1.45* 0.0212 * c^3 * 6e-4; % torsion constant of beam section [m^4], assume the thickness is t = 6mm
    GK = G*K;
    r = 120/2000;
    thetap = diff(theta)*pi/180;
    Mt = GK * thetap;
    Tx = Mt/(r);
    % A = 2*80*110*10^(-6);
    A = 2 * 54000*0.6*10^(-6);
    J = pi*(110^4 - 80^4)/32 * 10^(-12)
    % J = 2 * Ixx
    tau = Mt/J;

    Ultimate_norm = 633E6 ;
    Ultimate_shear = 31E6 ;
    % Ultimate_shear = 2015/J ;

    SF_n = Ultimate_norm / max(abs(sigma));
    SF_s = Ultimate_shear / max(abs(tau));




end