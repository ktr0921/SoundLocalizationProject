% NEW_FEEDBACK Initializes neural network for feedback LQR controller.
%
% Inputs:
%     kneuron_prev - Previous K
%     inputRate    - (V_IA) Auditory neuron output spike rate (before being weighted by J)
%
% Outputs:
%     kneuron - contains the output u and all the recursive information.
%
% K = diag[(1-2*eta*alpha)/2*eta*V_IA_1  (1-2*eta*alpha)/2*eta*V_IA_2 ... (1-2*eta*alpha)/2*eta*V_IA_n]
function kneuron = update_feedback(inputRate, kneuron_prev)
   % For the rate based neuron. The output rate is a function of the
   % excitatory input rate, i.e. G*w*R_i, where gain is defined by the 
   % inhibitory input rate V_IA. R_i is either the alpha as an input rate,
   % or the baseline rate (spontaneous activity) of the same neuron. 
   kneuron.w_input = kneuron_prev.w_input;%kneuron_prev.w_input;
   kneuron.baseline = kneuron_prev.baseline;   
%    kneuron.K = baseline./(w*VI_A);
   kneuron.K = kneuron.baseline./(kneuron.w_input * inputRate);
   
   if false % If spiking neuron
         r_e = kneuron_prev.baseline; % Maybe it can be an excitatory input rate
         r_i = inputRate% inhibitory input rate (V_IA)
         V = kneuron_prev.V; % Membrane potential

         shunting = g_e*tau_e*r_e + g_i*tau_i*r_i; % Average amount of shunting for excitatory and inhibitory inputs
         I_variance = g_E^2*tau_e*r_e*(V - Ee)^2 + g_i^2*tau_i*r_i*(V - E_i)^2
         I = kneuron_prev + I_variance;
         
         % This comes from Chance 2002. Need to pass it through a LIF
         % neuron to generate the spikes.
   end
end