% u = controller(counter, con_type, recon, jeff, Alpha, samp_per)
function control = controller(tstep, control, control_type, recon, supervisor, jeff, Alpha, dt, DEBUG)
   if strcmp('non', control_type)
      u = recon.V; % reconstruction neuron output is error btwn vis & aud
   elseif strcmp('lqr_ff', control_type)
      % u = K*(x - x_d) + u_d;
      V_IA = jeff.V_IA;
      r = supervisor.V_V;
      x_d = zeros(size(V_IA));
      u_d = zeros(size(V_IA));
      k = zeros(length(V_IA),1);
      
      qi = 1;
      Ri = 200;
      
      for i = 1:numel(jeff.V_IA)
          x_d(i) = (r(i)/sum(V_IA.^2))*V_IA(i); % x_d = inv(C)*r
          u_d(i) = -(r(i)/sum(V_IA.^2))*Alpha(i);% u_d = -inv(B)*A*x_d
          k(i) = (Alpha(i)/V_IA(i)) - qi*V_IA(i)/(Ri*2*Alpha(i)); % 
          k(i) = min(k(i),1e10);
      end
      K = diag(k);
      
%       inputRate = V_IA;
%       control.kn = update_feedback(inputRate, control.kn);
%       K = control.kn.K;
      
      for ni=1:numel(jeff.V_IA)
          x = jeff.J(:,ni);
          % u should be a 1x20 vector
          u(:,ni) = (K*(x_d - x) + u_d)'; 
      end
      
      
   elseif strcmp('lqr', control_type)
      [NA,NR] = size(jeff.J);
      eta = 1/jeff.timeconst;
      D     = 0;
      Q     = 100*ones(NA+1,NA+1);
      R     = 1;
      B     = -2*eta*jeff.V_IA;
      C     = B.';
      B_aug = [0 B];
      C_aug = B_aug.';        
      eta   = 1/jeff.timeconst;

      u     = zeros(1,NR);
      xI    = control.xI;
      for loop = 1:NR
        A      = -2*eta*diag(Alpha(:,loop));
%         A_aug  = [0 -C; zeros(1,NA) A];
        A_aug  = [0 zeros(1,NA); -C A];
        x      = jeff.J(:,loop);
        if all(abs(B)>eps)
           fprint('\n');
        end

%         [A_, B_, C_, T, K] = ctrbf(A_aug,B_aug',C_aug');
        %returns a decomposition into the controllable/uncontrollable subspaces.

        % The number of controllable states is sum(K)
%         Nc = sum(K);

        % The number of non-controllable states is then simply the
        % number of total states subtracted by the number of
        % controllable states.
%         Nnc = size(A_aug,1)-Nc;

        % Therefore, we need to inspect the A_{nc}, the matrix
        % responsible for the stability of the uncontrollable states.
        % This is the first sum(K) x sum(K) element of the matrix ABAR,
        % which is already in Kalman Decomposition Form.
%         A_nc = A_(1:Nnc,1:Nnc);

        % If A_nc is unstable, then the system in not stabilizable
   %         fprintf('%f\n',max(eig(A_nc)));
%         if max(eig(A_nc))>= 0
%             u(loop) = recon.V(loop);
% 
%         % Otherwise, it is stabilizable and so we can apply LQR
%         else
   %             display('Stabilizable system')
            K_lqr = lqr(A_aug,B_aug',Q,R);
%             if isempty(K_lqr)
%                K_lqr = zeros(1,26);
%             elseif all(K_lqr)==0
% %                fprintf('K is zero\n');
%             else               
%                ['K_lqr=[' fprintf('%.2g\t', K_lqr) ']' fprintf('\n')];
%             end
            %             xI(loop)= xI(loop) + recon.V(loop)*samp_per;
            xI(loop) = recon.V(loop);

            x_aug    = [xI(loop);x];
            u(loop)  = -K_lqr*x_aug;
%         end
      end
      control.xI = xI;
      
   else
     u = recon.V;
   end
   
   control.u  = u * 1e5;
%    if DEBUG
%        
%       control.DEBUG.u(DEBUG,:) = u;
%    end
end

%% Controller figures
% % Code to see the impact of K compared to Franosch error in recon.V
% % Alpha<=1e-3 is too small for weights to influence control
% Alpha  = 1e-2*ones(NA,NR);
% A      = -2*eta*diag(Alpha(:,loop));
% A_aug  = [0 zeros(1,NA); -C A];
% 
% Qvec = [1:10:101];
% Rvec = Qvec; Rvec(Rvec==0) = 1;
% nQ   = length(Qvec); 
% nR   = length(Rvec);
% c    = zeros(NR,nQ,nR);
% k    = zeros(NA+1,NR,nQ,nR);
% for qi=1:nQ
%    for ri=1:nR
%       Q=Qvec(qi)*ones(NA+1,NA+1); 
%       R=Rvec(ri); 
%       for ni=1:NR
%          k(:,ni,qi,ri) = lqr(A_aug,B_aug',Q,R); 
%          x_aug = [recon.V(ni); x]; 
%          c(ni,qi,ri) = -k(:,ni,qi,ri)'*x_aug; 
%       end
%    end
% end
% figure; 
% plot(reshape(c,[NR nQ*nR]));
% hold on; 
% plot(recon.V,'k');
% 
% figure; nr = ceil(sqrt(nR)); nc = ceil(nR/nr);
% for ri=1:nR
%    subplot(nr,nc,ri);
%    plot(c(:,:,ri));
%    hold on; 
%    plot(recon.V,'k');
%    title(sprintf('R=%d',Rvec(ri)));
%    ylabel('Control');
%    xlabel('Recon neuron');
% end
% legend(str2legend('Q=',Qvec));
% 
% figure; nr = ceil(sqrt(nQ)); nc = ceil(nQ/nr);
% for qi=1:nQ
%    subplot(nr,nc,qi);
%    plot(permute(c(:,qi,:),[1 3 2]));
%    hold on; 
%    plot(recon.V,'k');
%    title(sprintf('Q=%d',Qvec(qi)));
%    ylabel('Control');
%    xlabel('Recon neuron');
% end
% legend(str2legend('R=',Rvec));


