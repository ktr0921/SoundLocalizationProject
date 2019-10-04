% [v, spikes] = updateLIFNeuron(v, spikes, threshold, tau, t_refrac, v_reset, R, I, dt)
%
% Udpate LIF neuron potential, & determine if spike occurs or not
% Allows a vector of neurons to be updated simultaneously
% Inputs:
%  v        - matrix of neuron potentials from previous timestep (matrix)
%  spikes   - matrix of spiking history to implement refractory period (matrix)
%  tau      - decay time constant for neuron potential (scalar)
%  t_refrac - refractory period of neurons (scalar)
%  R        - membrane resistant (scalar)
%  I        - input to neurons (vector of size matching v, or scalar)
%  dt       - time resolution (scalar)
function [v, spikes] = updateLIFNeuron(v, spikes, threshold, tau, t_refrac, v_reset, R, I, dt)
     dv = (-v/tau + R*I)*dt;
     % Calculate how many timesteps back the refractory period lasts
     v  = v + dv;
     
     % Check number of timesteps of spike memory matches refractory period
     if size(spikes,1) > ceil(t_refrac/dt)
        spikes = spikes(end-ceil(t_refrac/dt)+1:end,:);
     elseif size(spikes,1) < ceil(t_refrac/dt - dt)
        fprintf('\twarning:updateLIFNeuron: insufficient spike memory to implement refractory period\n\n');
     end
     % number of spikes remembered is determined by refractory period
     % --> any spikes in a neuron's vector of spikes means it can't fire
     v  = (v > 0 & ~any(spikes)) .* v; % lower bound with 0
    
    % If the potential is at threshold level or higher, generate spike
    spikes = (v >= threshold);
    
    if any(spikes(:,end))
%        disp('pause');
    end
    
    % Reset the potential to v_r if the potential is above threshold.
    % Otherwise, don't do anything
    v = v.*(v < threshold) + v_reset * (v >= threshold); % NR x 1
    
end