function jeff = update_weight(tstep, jeff, mode, model, source, controller, dt, Alpha,DEBUG)
% The function jeff = update_weight(jeff, mode) allows the user to change 
% the weight (J) of the jeffress, by taking into account the learning rate 
% (tau), the input rate (before taking dot product with J, i.e. V_IA) and 
% the final output rate (after taking dor product with J, i.e. V_A). 
% There are two modes offerred in this function, unsupervised learning and
% supervised learning. Please input 'unsup' in mode for unsupervised 
% learning and 'sup' for supervised learning.
% If the mode is unsupervised, then there are two models which can be 
% controlled by the variable mode. The variable model is a string variable. 
% Input 'no' for a formula without offset (has been reduced by the base 
% rate and amplitude of the source) and input 'o' for a formula with offset
% If the mode is supervised, then insert 'stdp' for the model variable to
% apply the equation from Franosch 2013 MIT paper.
% Alpha is for the component of A matrix of the system equation 
% Jdot = A*J + V_IA*w, where the equation is applied once per column (split
% the column and then apply the equation one by one)
      u = controller.u; % When not controlled, this is the Visual Spike Rate
    tau = jeff.timeconst;
     NA = length(jeff.neuron);
     NR = size(jeff.J,2);
    vJ = jeff.V_IA;%     vJ  = jeff.V_IA;   % Jeffress ladder output values (vector, 1 x NJ)
    vR = u;%      vR  = jeff.V_A;    % Reconstruction neuron output values (vector, NR x 1)
    if mode == 'uns'
%         if model == 'no'
%             base = source.baseline;
%             ampl = source.amplitude;
%             vi = vi - (base + ampl);
%             vj = vj - (base + ampl);
%         end
%         
        % Now there is a value of J for each relationship between pre-synaptic
        % neuron (NJ of them) and post-synaptic neuron (NR of them). Therefore,
        % there are NJ*NR change in total. Look at the book for how it boils
        % down to just a simple matrix multiplication.
      dJ        = -1/jeff.timeconst * vJ*u;
      jeff.J    = jeff.J + dJ*dt;
      jeff.Jdot = dJ;
      
    elseif mode == 'sup'
        if model == 'stdp_var'
            % Update weights for a vector of auditory neurons to a
            % single reconstruction neuron - ie loop through reconstrution neurons
            A = diag(Alpha(1,:));
            for ni=1:NA
%                 x    = jeff.J(ni,:);
%                 Jdot = -1/jeff.timeconst * ((x*A) + vJ(ni)*u); % control vector
%                 x    = jeff.J(ni,:) + Jdot*dt;
%                 x(x>jeff.max_wgt) = jeff.max_wgt;
%                 x(x<jeff.min_wgt) = jeff.min_wgt;
%                 jeff.J(ni,:)      = x;
%                 jeff.Jdot(ni,:)   = Jdot;

                if isfield(controller,'u_d'), u = controller.u(:,ni)';end % If using LQR ff-fb controller, u is not a vector, but a matrix
                x    = jeff.J(:,ni);
                Jdot = -1/jeff.timeconst * ((A*x) + vJ(ni)*u'); % control vector
                x    = jeff.J(:,ni) + Jdot*dt;
                x(x>jeff.max_wgt) = jeff.max_wgt;
                x(x<jeff.min_wgt) = jeff.min_wgt;
                jeff.J(:,ni)      = x;
                jeff.Jdot(:,ni)   = Jdot;
            end
            if DEBUG
               % If debugging, log Jdot contributors
                jeff.DEBUG.B(DEBUG,:)        = jeff.V_IA;
                jeff.DEBUG.u(DEBUG,:)        = u;
            end
            
        elseif model == 'stdp_fra'
            % Apply Franosch STDP Formula
            gamma     = 1e-3;
            % it may take ages for correlation estimate to converge to
            % correct scale, so initialise it with first estimation
            if ~any(jeff.corr(:))
               jeff.corr = jeff.V_IA'*u;
            end
            jeff.corr = jeff.corr + jeff.V_IA'*u*dt;
            dJ        = -2/jeff.timeconst * (gamma*jeff.J + jeff.corr);
            x         = jeff.J + dJ*dt;
            x(x>jeff.max_wgt) = jeff.max_wgt;
            x(x<jeff.min_wgt) = jeff.min_wgt;
            jeff.J    = x;
            jeff.Jdot = dJ;

            if DEBUG
               jeff.DEBUG.corr(DEBUG,:,:) = jeff.corr;
            end
        end
    end
    
   if DEBUG
      % If debugging, log Jdot contributors
       jeff.DEBUG.Jdot(DEBUG,:,:) = jeff.Jdot;
       jeff.DEBUG.J(DEBUG,:,:)    = jeff.J;
   end
end