% jeff = update_jeffress(tstep, dt, jeff, source, tt, epsp, window_sample, basemodel, DEBUG)
function jeff = update_jeffress(tstep, dt, jeff, source, epsp, neuronType, DEBUG,base)
   % DEBUG is either false, or the debug sample number
   %
   % Update jeffress ladder neuron potentials using the threshold detector (V_IA)
   % of the Jeffress Ladder 'jeff', with the auditory source potential
   % coming to the Jeffress Ladder from the left & right ears.
   % Update the Jeffress Ladder at time tt.
   tt = tstep*dt;
   threshold = jeff.threshold;
   NJ = length(jeff.neuron);

   % if neurontype is rate based then propagate sin wave down the ladder,
   % if it's spike based, propagate a spike down the ladder
   if strcmp('r',neuronType)
       V_left  = zeros(1,NJ);
       V_right = zeros(1,NJ);
   elseif strcmp('s',neuronType)
       % push old samples down for windowed data
       [jeff.Lspikes, jeff.Rspikes, jeff.spikes] = ageData(jeff.Lspikes, jeff.Rspikes, jeff.spikes);
       V_left  = 0;
       V_right = 0;
   end

   for nsource = 1:length(source)
       t_src   = source(nsource).t_source;
       sigma   = source(nsource).gauss_std;

       theta_s  = source(nsource).location;
       A        = source(nsource).amplitude;
       f        = source(nsource).frequency;
       omega    = 2*pi*f;
       baseline = base;%source(nsource).baseline;

       source_delay_L = (-theta_s + pi/2)/pi*jeff.delta_lr;    % delay to ladder - scalar (see derivation in book)
       source_delay_R = (theta_s + pi/2)/pi*jeff.delta_lr;     % delay to ladder - scalar (see derivation in book)
       propagation_delay_L = (jeff.neuron + pi/2)/pi*jeff.delta_lr;    % delay to neuron in ladder - NJ x 1
       propagation_delay_R = (pi/2 - jeff.neuron)/pi*jeff.delta_lr;	% delay to neuron in ladder - NJ x 1

       if neuronType=='r'
           delay_L = source_delay_L + propagation_delay_L;	% total delay - NJ x 1
           delay_R = source_delay_R + propagation_delay_R; % total delay - NJ x 1

           % Don't rectify left & right separately, only the sum
           V_left  = V_left  + (A*cos(omega*(tt-delay_L))').*heaviside(tt-delay_L)' + baseline + normrnd(0,sigma);
           V_right = V_right + (A*cos(omega*(tt-delay_R))').*heaviside(tt-delay_R)' + baseline + normrnd(0,sigma);
%            V_left  = V_left  + (A*sin(omega*(tt-delay_L))') + baseline + normrnd(0,sigma);
%            V_right = V_right + (A*sin(omega*(tt-delay_R))') + baseline + normrnd(0,sigma);

           if DEBUG
               delay_L = source_delay_L + (0:(NJ-1))'/(NJ-1)*jeff.delta_lr;
               delay_R = source_delay_R + ((NJ-1):-1:0)'/(NJ-1)*jeff.delta_lr;
           end


       elseif neuronType=='s'
           % For spiking input, spike on left & right when input arrives, &
           % then propagate the spike down the ladder. Total input into left
           % & right is addition of all sources, & then determine spiking

           % ITD for signals < 2KHz is: tau = (a/c)2sin(theta)
           % ITD for signals > 2KHz is: tau = (a/c)(theta + sin(theta))
           % a = head radius, c = speed of sound, theta wrt 0 at left ear
           a   = 85.5e-3; c = 340*16; % a [m], c [m/s]
           ITD = ternaryOp(f<=2e3, (a/c)*2*sin(theta_s), (a/c)*((theta_s) + sin(theta_s)) );
           % ITD<0 for theta_s<0, ITD>0 for theta_s>0
           if ITD<0
               source_delay_L = 0; source_delay_R = -ITD;
           elseif ITD>0
               source_delay_R = 0; source_delay_L = +ITD;
           elseif ITD==0
               source_delay_R = 0; source_delay_L = 0;
           end
           V_left   = V_left  + heaviside((tt-t_src)-source_delay_L).*A*(cos(omega*((tt-t_src)-source_delay_L)))' + baseline + normrnd(0,sigma);
           V_right  = V_right + heaviside((tt-t_src)-source_delay_R).*A*(cos(omega*((tt-t_src)-source_delay_R)))' + baseline + normrnd(0,sigma);
       end
   end % end for each source

   if neuronType=='r' % rate based neurons in ladder
       % Full wave rectification
       jeff.V_L = bsxfun(@max, V_left, 0);  % add new sample - NJ x 1 (Jeffress potential before summing)
       jeff.V_R = bsxfun(@max, V_right, 0); % NJ x 1 (Jeffress potential before summing)
%        jeff.V_L = V_left;
%        jeff.V_R = V_right;
  
       I = jeff.V_L' + jeff.V_R'; % NJ x 1: Jeffress input potential, skipping the epsp
%        I = conv(I,epsp,'same');
       I(I < jeff.threshold + (baseline * 2)) = baseline; % I = (I > jeff.threshold).*I; % apply if larger than threshold, otherwise it is 0.
       jeff.V_IA = I;

       % Update final auditory neuron output, that is input to recon neurons
       % (weighted spike rate [V_IA] is the input [V_A] to the reconstruction)
       
       jeff.V_A  = (jeff.J * I)'; % NR x 1 (contribution from Jeffress potential in reconstruction neuron)

   elseif neuronType=='s' % spiking neurons in ladder
       tau_m = jeff.tau_m; t_refrac = jeff.t_refrac; R_LIF = jeff.R_LIF; v_reset = jeff.v_reset;
       [jeff.V_L, jeff.Lspikes(end,:)] = updateLIFNeuron(jeff.V_L, jeff.Lspikes, threshold, tau_m, t_refrac, v_reset, R_LIF, V_left, dt);
       [jeff.V_R, jeff.Rspikes(end,:)] = updateLIFNeuron(jeff.V_R, jeff.Rspikes, threshold, tau_m, t_refrac, v_reset, R_LIF, V_right, dt);
       % need to record spike time so we can delay its propagation down ladder
       if any(jeff.Lspikes(end,:))
           jeff.Llastspike = tt;
       end
       if any(jeff.Rspikes(end,:))
           jeff.Rlastspike = tt;
       end
       delay_L  = jeff.Llastspike + propagation_delay_L;
       delay_R  = jeff.Rlastspike + propagation_delay_R;
       sympref('HeavisideAtOrigin', 1); % make heaviside 0 for input 0
       V_LA      = (heaviside((tt+dt-delay_L)) .* heaviside((delay_L+dt-tt)) .* any(jeff.Lspikes) )';
       V_RA      = (heaviside((tt+dt-delay_R)) .* heaviside((delay_R+dt-tt)) .* any(jeff.Rspikes) )';
       if any(V_LA) && any(V_RA)
       end
       I         = V_LA + V_RA;

       % Ladder neurons need both left & right spikes to fire --> threshold is twice as high
       [jeff.V_LIF, jeff.spikes(end,:)] = updateLIFNeuron(jeff.V_LIF, jeff.spikes, 0.002, tau_m, t_refrac, v_reset, R_LIF, I, dt);
       %       [jeff.V_LIF, jeff.spikes(end,:)] = updateLIFNeuron(jeff.V_LIF, jeff.spikes, 2*threshold, tau_m, t_refrac, v_reset, R_LIF, I, dt);
       jeff.V_IA = jeff.spikes(end,:); % input in to jeffress ladder neurons

       % Now convolve once again with the epsp to finally obtain V_A - this
       % takes output of jeffress neuron & converts to input of reconstruction
       % neuron.
       % Matrix of weights J * vector of jeffress neurons = vector of audio
       jeff.V_A = (jeff.J'*(jeff.spikes'*(flipud(epsp)*dt)))';


   end

   % if DEBUG && tt>0 && (round(tt/1e-4) == (tt/1e-4))
   %    figure;
   %    subplot(2,2,1), plot(jeff.DEBUG.VL([1 NJ],:)');
   %    subplot(2,2,2), plot(jeff.DEBUG.VR([1 NJ],:)');
   %    subplot(2,2,3), plot(jeff.DEBUG.VL([1 NJ],:)'+jeff.DEBUG.VR([1 NJ],:)');
   %    subplot(2,2,4), plot(V_left); hold on; plot(V_right); plot(conv_L); plot(conv_R); plot(conv_L + conv_R);
   %    lh=legend('L','R','L*\epsilon','R*\epsilon','(L+R)*\epsilon','topwest');
   % end

   if DEBUG
       jeff.DEBUG.V_A(DEBUG,:)    = jeff.V_A;
       jeff.DEBUG.V_IA(DEBUG,:)   = jeff.V_IA;
       jeff.DEBUG.V_L(DEBUG,:)    = jeff.V_L(end,:);
       jeff.DEBUG.V_R(DEBUG,:)    = jeff.V_R(end,:);
       if neuronType=='s'
           jeff.DEBUG.V_left(DEBUG,:) = V_left;
           jeff.DEBUG.V_right(DEBUG,:)= V_right;

           jeff.DEBUG.V_LIF(DEBUG,:)  = jeff.V_LIF;
           jeff.DEBUG.spikes(DEBUG,:) = jeff.spikes(end,:);
           jeff.DEBUG.Lspikes(DEBUG,:)= jeff.spikes(end,:);
           jeff.DEBUG.Rspikes(DEBUG,:)= jeff.spikes(end,:);
       end
   end % end debug
end % end function

% Push data down a timestep so current data sample can evolve
% Can operate on any number of inputs
function varargout = ageData(varargin)
for ni=1:length(varargin)
    varargin{ni}(1:end-1,:) = varargin{ni}(2:end,:);
end
varargout = varargin;
end
