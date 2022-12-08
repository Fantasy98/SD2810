function [ucrit,pcrit,zcrit,pconv,uvec] = flutter(M,K,Qip,neig,iter = 100)
% [ucrit,pcrit,zcrit] = flutter(M,K,Qip,neig)
% pconv : All laplacian frequency computed
% Compute flutter speed, frequency and critical mode.
%
% M : mass matrix
% K : stiffness matrix
% Qip : aerodynamic loads interpolation struct
% neig: number of eigenvector you want
% (c) 2004-2016 Dan Borglund <dodde@kth.se> and David Eller <dlr@kth.se>
  u = 14; 
  i = 0.1;
  % Since we know there are only 2 
  
  for iu = 1:iter
      
      [kbounds] = pk_bounds(u,M,K,Qip,neig);    
      for ieig = 1:neig
          [pu,vu] = pk_bisect(ieig,u,M,K,Qip,kbounds);
          pconv(ieig,iu) = pu;
          
          Vu(ieig,iu,:) = vu; 
            
      end
      
      u = u + i;
      uvec(iu) = u;
  end
  
  ucrit = max(uvec);
  pcrit = 0;
  zcrit = 0;
  % Sort out all real part >0
  indx = real(pconv) > 0;
  for imode = 1:neig
      
      % fprintf("At mode %.2d \n",imode);
      if length(unique(indx(imode,:))) > 1
          % fprintf("Found flutter at mode %.2d\n",imode)
          loc = min( find(indx(imode,:)==1) );
          % loc = loc_list(1);    
          u_choose = uvec(loc);
          % fprintf("Found flutter Speed %.2d\n",u_choose);
          if ucrit >= u_choose
            u_save = u_choose
            ucrit = u_choose; 
            p_flutter = pconv(imode,loc);
            pcrit = imag(p_flutter) * ucrit/ (2*pi*Qip.bref);
            zcrit = Vu(imode,loc,:);
          else
            ucrit = u_save;
            pcrit = 0;
            zcrit = 0;
           end
 
  end





  % % %   % Find if the real part of phat >0
  %   loc_list = find( real(pconv(imode,:))>0 );
  % % %   % If exist, the number of coloum will no longer be 0
  %    if size(loc_list,2) > 0
  %       fprintf("Find the fultter mode is at mode  %.d \n",imode);
  %       nmode = imode;
        
  % %       %the first loctation is flutter speed
  %       loc = loc_list(1);
  % %       % p_flutter = pconv(imode,loc);
  % %       % loc1 = 0;
            
  %   end
  % end
  % ucrit = 10;
  % pcrit = 10;
  % zcrit = 10;
  % set return values
  % ucrit = uvec(loc);
  % pcrit = imag(p_flutter) * ucrit/ (2*pi*Qip.bref);
  % zcrit = squeeze(Vu(nmode,loc,:));

end
