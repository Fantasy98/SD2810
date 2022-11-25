function [ucrit,pcrit,zcrit] = flutter(M,K,Qip)
% [ucrit,pcrit,zcrit] = flutter(M,K,Qip)
%
% Compute flutter speed, frequency and critical mode.
%
% M : mass matrix
% K : stiffness matrix
% Qip : aerodynamic loads interpolation struct
%
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>
  u = 15.5; 
  i = 0.05;
  % Since we know there are only 2 
  neig = 2;
  for iu = 1:100
      
      [kbounds] = pk_bounds(u,M,K,Qip,neig);    
      for ieig = 1:neig
          [pu,vu] = pk_bisect(ieig,u,M,K,Qip,kbounds);
          pconv(ieig,iu) = pu;
          Vu(ieig,iu,:) = vu;
      end
      
      u = u + i;
      uvec(iu) = u;
  end
  for imode = 1:neig
    % Find if the real part of phat >0
    loc_list = find( real(pconv(imode,:))>0 );
    % If exist, the number of coloum will no longer be 0
    if size(loc_list,2) > 0
        fprintf("Find the fultter mode is at mode  %.d \n",imode);
        nmode = imode;
        % Corresponding the first loctation is flutter speed
        loc = loc_list(1);
        p_flutter = pconv(imode,loc)
        % loc1 = 0;
            
    end
  end

  % set return values
  ucrit = uvec(loc);
  pcrit = imag(p_flutter) * ucrit/ (2*pi*Qip.bref);
  zcrit = Vu(nmode,loc,:);

end
