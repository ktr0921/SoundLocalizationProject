function visual = new_supervisor(NV,type, gaussian_data, sigma_F, recon, neuronType, window_size, DEBUG, NOUT)
% The function supervisor = pl_new_supervisor(n_supervisor,type,K) creates
% a new supervisor neural systems with n_supervisor neurons.
% As these neurons are connected to the reconstruction neurons, the
% variable Fi represents the weight given between the visual neuron and the
% reconstruction neuron. For example, Fi(1,4) is the weight given from neuron
% 1 of the supervisor sensor and the fourth neuron of the reconstruction
% neuron. See Franosch et al., 2013 for this model.
% The type represents the spacing between the neurons in the supervisor
% sensor (i.e. the visual sensor). Make 'u' as the input if uniform spacing
% is desired. Otherwise, type in 'g' for gaussian spacing.
% When Gaussian type is chosen, the user must give the guassian_data input.
% Make the first entry as the mean of the Gaussian and the second as the
% standard deviation of the Gaussian distribution.
          % Initialise the spikes xtilde

   NR = length(recon.neuron);
   Fi = zeros(NV,NR);
   if type == 'l'
       visual.supervisor = linspace(-pi/2, pi/2, NV);
       visual.supervisor = visual.supervisor';
       for nid = 1:NV
           theta_i = visual.supervisor(nid); % angle of individ visual neuron  NK x 1
           Theta_j = recon.neuron.';         % angle vector of reconstruction neurons NR x 1
           if sigma_F==0
               % If std dev is 0 set closest visual neuron to have rate of A (amp),
               % and all other visual neurons to have a rate of 0
               d = Theta_j - theta_i; % get diff btwn recon & visual angles
               [~,i] = min(abs(d)); % index of visual neuron closest to recon neuron
               Fi(nid,:) = zeros(size(Theta_j));
               Fi(nid,i) = 1;
           else
               scaling(nid) = normpdf(theta_i, theta_i, sigma_F);
               Fi(nid,:) = normpdf(Theta_j,theta_i,sigma_F)/scaling(nid);
           end
       end
       % make sure all reconstruction neurons are being used
       if any(all(Fi==0))
         unused = find(all(Fi==0));
         copy   = unused-1; 
         copy(copy<1) = NR + copy(copy<1); % circular copy in case 1st neuron unused
         % copy recon neuron before the one that hasn't been used
         Fi(:,unused) = Fi(:,copy);
       end
       % scale each column by the same amount to not artificially increase
       % edge values compared to more central values
       Fi = Fi./repmat(max(sum(Fi)),[NV 1]); % ensure wgts sum to 1 so not adding power

       % Determine F, the weight between the visual and reconstruction neurons
       % F_ij is the weight between the i^th visual neuron and the j^th reconstruction neuron
       % z.F = eye(NK,NR); 
       visual.F = Fi;
   elseif type == 'g'
       mu = gaussian_data(1);
       sigma = gaussian_data(2);
       start = normcdf(-pi/2,mu,sigma);
       finish = normcdf(pi/2,mu,sigma);
       increment = (finish-start)/(n_auditory-1);
       currentprob = start;
       for i = 1:NJ
           visual.supervisor(i) = norminv(currentprob, mu, sigma);
           currentprob = currentprob + increment;
       end
       if visual.supervisor(end)~=pi/2
           visual.supervisor(end) = pi/2;
       end
       visual.supervisor = visual.supervisor';
       % z.F = eye(NK,NR); 
       for nid = 1:NV
           theta_i = visual.supervisor(nid);
           theta_j = recon.neuron.';
           scaling = normpdf(theta_i, theta_i, sigma_F);
           Fi(nid,:) = normpdf(theta_j,theta_i,sigma_F)/scaling;
       end
       visual.F = Fi;
   else
       display('type must be either ''l'' or ''g''');
       return
   end

   if neuronType=='s'
      visual.V_LIF  = zeros(1, NV);
      visual.spikes = zeros(window_size, NV); % window x NV
   end

   if DEBUG
      visual.DEBUG.V_LIF  = zeros(NOUT, NV); 
      visual.DEBUG.spikes = zeros(NOUT, NV); 
      visual.DEBUG.V_IV   = zeros(NOUT, NV); 
      visual.DEBUG.V_V    = zeros(NOUT, NR); 
   end
end