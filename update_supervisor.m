% visual = update_supervisor(counter,visual,source,samp_per,epsp,neuronType,DEBUG)
% Updates supervisor rate, V_V - uses Poisson or LIF neurons
function visual = update_supervisor(tstep, dt, visual, source, epsp, neuronType, DEBUG)
   NK = length(visual.supervisor);
   NV = length(visual.supervisor);
   rate_visual = zeros(1, NV);
   for nsource = 1:length(source)
       % extract source location and convert such that the angle is from -pi
       % to pi
       A = source(nsource).visualAmp;
       theta = source(nsource).location;
       theta_std = source(nsource).loc_std;
       %     rate_source = source.rate;
       while theta<-pi
           theta = theta+2*pi;
       end
       while theta>pi
           theta = theta-2*pi;
       end

       % Apply the model derived from the notebook
       neurons = visual.supervisor;
       if theta_std==0
           % If std dev is 0 set closest visual neuron to have rate of A (amp),
           % and all other visual neurons to have a rate of 0
           d = neurons - theta; % difference btwn visual neuron angle & source angle
           [~,i] = min(abs(d)); % index of visual neuron closest to source angle
   %         rate_visual = zeros(size(neurons)); % set all visual neuron rates to 0
           rate_visual(i) = A;  % set closest visual neuron closet to have rate of amps
       else
           f_theta = normpdf(neurons, theta, theta_std);
           rate_visual = rate_visual + A*f_theta/normpdf(theta, theta, theta_std);
       end
   end

   % TODO; move these parameters if we keep the spiking method
   visualNeuron = 'LIF';
   tau_v        = 1e-3; 
   threshold_V  = 1; 
   t_refrac     = dt;
   v_r          = 0;
   R_LIF        = 1e4; 
   v_reset      = 0;
   
   visual.V_IV = rate_visual;                      % source input
   if neuronType == 'r'                             % if it is rate based
%        visual.V_V = (visual.F).'*visual.V_IV;      % the potential is straight forward
       visual.V_V = visual.V_IV*visual.F;      % the potential is straight forward
       
   elseif neuronType =='s'                          % if it is spike-based
       % Spike according to Poisson probability 
       if strcmpi(visualNeuron,'Poisson')
          % Convert the rate from spikes/sec to spikes/sample
          rate_s = rate_visual*dt;
          % Now generate a random number from uniformly distributed random
          % variable. If that number is less than rate_s, we create a spike. The
          % ensemble average spike will then be rate_s at that particular time.
          uniform_rv = rand(size(rate_s));

          % Create a spike everytime the RV generated is less than rate_s.
          % The new data is concantenated to the right and the first element is
          % chopped off (visualise by using window)
          visual.spikes = [visual.spikes (uniform_rv <= rate_s)*1.0/dt];   % update the most recent information
          visual.spikes(:,1) = [];                   % chop off the first element

      elseif strcmpi(visualNeuron,'LIF')
         visual.spikes = ageData(visual.spikes);
         % Apply the dynamical equation of LIF
         [visual.V_LIF, visual.spikes(end,:)] = updateLIFNeuron(visual.V_LIF, visual.spikes, threshold_V, tau_v, 0, v_reset, R_LIF, rate_visual, dt);
         if DEBUG
            visual.DEBUG.V_LIF(DEBUG,:)  = visual.V_LIF(end,:); 
            visual.DEBUG.spikes(DEBUG,:) = visual.spikes(end,:); 
         end
       end
       % Now convolute it with the epsp and weigh by F to get V_V
       visual.V_V = (visual.F' * visual.spikes'*(flipud(epsp))*dt)';
   end

   if DEBUG
      visual.DEBUG.V_IV(DEBUG,:) = visual.V_IV;
      visual.DEBUG.V_V(DEBUG,:)  = visual.V_V;
   end

end

% Push data down a timestep so current data sample can evolve
% Can operate on any number of inputs 
function varargout = ageData(varargin)
   for ni=1:length(varargin)
      varargin{ni}(1:end-1,:) = varargin{ni}(2:end,:);
   end
   varargout = varargin;
end