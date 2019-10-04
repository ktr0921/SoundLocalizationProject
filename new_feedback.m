% NEW_FEEDBACK Initializes neural network for feedback LQR controller.
%
% Inputs:
%     K     - First K
%     eta   - (scalar)
%     alpha - (vector)
%
% Outputs:
%     controller - contains the output u and all the recursive information.
%
% K = diag[(1-2*eta*alpha)/2*eta*V_IA_1  (1-2*eta*alpha)/2*eta*V_IA_2 ... (1-2*eta*alpha)/2*eta*V_IA_n]
function kneuron = new_feedback(NA, eta, alpha)
   kneuron.baseline = 1-2*eta*alpha; % Spike rate of spontaneous activity
   kneuron.K = kneuron.baseline * diag(0.1*ones(NA,1)); % K with no input. Just the baseline rate. K  is the output spike rate of the neuron
   kneuron.w_input = ones(1,NA)*2*eta; % Check the dimension conicides with V_IA.
end