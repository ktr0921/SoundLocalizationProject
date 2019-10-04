function jeff = new_jeffress(NA, NR, type, model, gaussian_data, delta_lri, tau, neuronType, max_wgt, min_wgt, threshold, window, R_LIF, tau_m, v_reset, t_refrac, DEBUG, NOUT)
%  jeff = new_jeffress(NJ, type, gaussian_data, J_ij, delta_lri, tau, ...
%                      neuronType, max_wgt, min_wgt, threshold, window, ...
%                      R_LIF, tau_m, v_reset, DEBUG, NOUT)
% receives the number of neurons in the Jeffress ladder (NJ) and the
% density type of the ladder. Its output is the neurons' angles involved in
% the Jeffress ladder and the time required for the signal to travel from
% left side of the ladder to the right side.
% For a linear density accross all angles from -pi/2 to pi/2 radians,
% type in 'l' (for linear) for the type. Leave gaussian_data as [] when
% using this
% For a Gaussian-like density along -pi/2 to pi/2 radians, type in a
% vector with 'g' as the type. The gaussian_data is mean as the first
% element, then its standard deviation as second
% delta_lri is the time taken for the signal to go across from the end of
% left ladder to the end of the right side of Jeffress Ladder.
% tau is the learning rate. In unsupervised learning, this is the learning
% rate (tau), while in supervised learning, this is taken to be eta, the
% learning constant (see stdp).
% max_wgt and min_wgt provide upper and lower bounds for plastic weights.
      
   jeff.J = rand(NA,NR)*0.2;%(1/NR); % wgt btwn i^th Jeffress neuron and the j^th recon neuron
   if type(1) == 'l'
       jeff.neuron    = linspace(-pi/2, pi/2, NA);
       jeff.neuron    = jeff.neuron';
       jeff.delta_lr  = delta_lri;
       jeff.timeconst = tau;
       jeff.threshold = threshold;
       
   elseif type(1) == 'g'
       jeff.threshold = threshold;
       jeff.delta_lr  = delta_lri;
       jeff.timeconst = tau;
       mu             = gaussian_data(1);
       sigma          = gaussian_data(2);
       start          = normcdf(-pi/2,mu,sigma);
       finish         = normcdf(pi/2,mu,sigma);
       increment      = (finish-start)/(NA-1);
       currentprob    = start;
       for i = 1:NA
           jeff.neuron(i) = norminv(currentprob, mu, sigma);
           currentprob = currentprob + increment;
       end
       if jeff.neuron(end)~=pi/2
           jeff.neuron(end) = pi/2;
       end
       jeff.neuron = jeff.neuron';
   else
       display('type must be either ''l'' or ''g''');
       return
   end
  
   if neuronType=='r' % rate based neuron
     jeff.V_L    = zeros(1, NA); % initialise left jeffress input
     jeff.V_R    = zeros(1, NA); % initialise right jeffress input
     jeff.maxI = [];
     
     % Following lines are to prevent an error on update_jeffres
     jeff.Lspikes = zeros(window, 1); 
     jeff.Rspikes = zeros(window, 1); 
     jeff.Llastspike  = Inf; % left spike time
     jeff.Rlastspike  = Inf; % left spike time
   elseif neuronType=='s' % spiking neuron
     jeff.V_L    = 0; % initialise left jeffress input
     jeff.V_R    = 0; % initialise left jeffress input
     % Update V_LIF according to the LIF dynamical system equation
     jeff.V_LIF  = zeros(1, NA);      % Initialise V_LIF
     jeff.spikes = zeros(window, NA); % NR x window sample, spike after thresholding
      
     jeff.Lspikes = zeros(window, 1); 
     jeff.Rspikes = zeros(window, 1); 
     jeff.Llastspike  = Inf; % left spike time
     jeff.Rlastspike  = Inf; % left spike time

     jeff.R_LIF   = R_LIF;
     jeff.tau_m   = tau_m; 
     jeff.v_reset = v_reset; 
     jeff.t_refrac= t_refrac;
   end

   %TODO: J_ij must be column vector i.e. NJ by 1
   jeff.max_wgt = max_wgt;
   jeff.min_wgt = min_wgt;
   
   if model=='stdp_fra'
      % In STDP, we are required to integrate. Thus, we need to
      % save the previous value
        jeff.corr = zeros(NA,NR);  % Initialise the correlation term
   end
    
   if DEBUG
     jeff.DEBUG.B    = zeros(NOUT, NA); 
     jeff.DEBUG.u    = zeros(NOUT, NR); 
     jeff.DEBUG.J    = zeros(NOUT, NA, NR); 
     jeff.DEBUG.Jdot = zeros(NOUT, NA, NR); 

     jeff.DEBUG.V_A    = zeros(NOUT, NR); 
     jeff.DEBUG.V_IA   = zeros(NOUT, NA); 
     if neuronType=='s'
        jeff.DEBUG.V_LIF  = zeros(NOUT, NA); 
        jeff.DEBUG.spikes = zeros(NOUT, NA);      
        jeff.DEBUG.Rspikes= zeros(NOUT, NA); 
        jeff.DEBUG.Lspikes= zeros(NOUT, NA); 
        jeff.DEBUG.V_L    = zeros(NOUT, 1); 
        jeff.DEBUG.V_R    = zeros(NOUT, 1); 
        jeff.DEBUG.V_left = zeros(NOUT, 1); 
        jeff.DEBUG.V_right= zeros(NOUT, 1); 
        
     elseif neuronType=='r'
        jeff.DEBUG.V_L    = zeros(NOUT, NA); 
        jeff.DEBUG.V_R    = zeros(NOUT, NA); 
     end
     if model=='stdp_fra'
        jeff.DEBUG.corr   = zeros(NOUT, NA, NR); 
     end
   end
end