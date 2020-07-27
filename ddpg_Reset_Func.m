function [InitialObservation, LoggedSignals] = ddpg_Reset_Func()
% Reset function to place custom cart-pole environment into a random
% initial state.

configuration; % Configuration
NOUT           = ceil(length(t)/out_dt);

% Create a new Jeffress Ladder
% Create reconstruction neurons
% New controller
% Initialise all variables which we want to monitor
jeff           = new_jeffress(NA, NR, 'l', 'stdp_fra', [], delta_lr, ...
    tau, neuronType, max_wgt, min_wgt, threshold, window_size, R_LIF, ...
    tau_m, v_reset, t_refrac, DEBUG, NOUT);
rec_Neuron     = new_recon(NR, 'l', [], DEBUG, NOUT);
visual         = new_supervisor(NV,'l',[],sigma_F, rec_Neuron, ...
    neuronType, window_size, DEBUG, NOUT);

location       = cell(length(t),1); % source location

% initialise source data. There are Nsource sound sources, so initialise 
% them one by one.
source(1:Nsource) = struct_source();
for nsource = 1:Nsource
    source(nsource).frequency = freq;
    source(nsource).amplitude = ampl_A;
    source(nsource).visualAmp = ampl_V;
    source(nsource).location  = lctn;
    source(nsource).loc_std   = mastersource.loc_std;
    source(nsource).baseline  = mastersource.baseline;
    source(nsource).gauss_std = mastersource.gauss_std;
    source(nsource).t_source  = mastersource.t_source;
    source(nsource).infchange = false;
end

% Create an array of index which helps us determine which data to use
% from the variable mastersource. The variable mastersource consists of
% a collection of source for all time, whereas the variable source
% consists of source data required for a particular time instance.
index = ones(Nsource,1);

tstep = 1;
if mod(tstep, out_dt)==0
    OUT = tstep/out_dt;
else
    OUT = false;
end

for nsource = 1:Nsource
    % Extract the correct source data from the mastersource.
    % Mastersource variable is an array of struct. There are
    % Nsource elements of mastersource, where each element consists
    % of a struct.
    source(nsource)  = struct_source(t(tstep), mastersource(nsource), index(nsource), source(nsource));
    % The infchange indicates that there is an information change
    % in the source. Therefore, this means that at this time, the
    % correct data that we need to extract from mastersource has
    % changed. Because it is already sorted such that the time is
    % increasing, then we can just simply increase the index for
    % that particular sound source.
    if source(nsource).infchange
        jeff.corr = zeros(NA,NR); % Initialise the correlation term
        index(nsource) = index(nsource) + 1;
    end
    location{tstep}(nsource) = source.location;
end

visual      = update_supervisor(tstep, dt, visual, source, epsp, neuronType, OUT);
jeff        = update_jeffress(tstep, dt, jeff, source, epsp, neuronType, OUT, 1e-100);
rec_Neuron  = update_reconstruction('sup', rec_Neuron, jeff, visual, inhibition, OUT);

LoggedSignals.jeff = jeff;
LoggedSignals.rec_Neuron = rec_Neuron;
LoggedSignals.source = source;
LoggedSignals.visual = visual;
LoggedSignals.location = location;
LoggedSignals.index = index;
LoggedSignals.tstep = tstep;
LoggedSignals.r = 0;
% LoggedSignals = [ddpg_Jeff(jeff, 'struct2mat'); rec_Neuron.neuron; source; visual; location; index];

% LoggedSignals.rec_Neuron.V = zeros(1, 20);
InitialObservation = LoggedSignals.rec_Neuron.V';

end