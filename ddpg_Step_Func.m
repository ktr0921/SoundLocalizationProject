function [NextObs, Reward, IsDone, LoggedSignals] = ddpg_Step_Func(Action, LoggedSignals)

configuration;
u_DDPG = struct();
u_DDPG.u = Action';

jeff       = LoggedSignals.jeff;
rec_Neuron = LoggedSignals.rec_Neuron;
source     = LoggedSignals.source;
visual     = LoggedSignals.visual;
location   = LoggedSignals.location;
index      = LoggedSignals.index;
tstep      = LoggedSignals.tstep;
r          = LoggedSignals.r;

if mod(tstep, out_dt)==0
    OUT = tstep/out_dt;
else
    OUT = false;
end

Alpha = 1e-4*ones(NA,NR);
jeff = update_weight(tstep, jeff, 'sup', 'stdp_var', source, u_DDPG, dt, Alpha, OUT);

tstep = tstep + 1;
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
LoggedSignals.r = r - sum(abs(rec_Neuron.V));

% LoggedSignals = [ddpg_Jeff(jeff, 'struct2mat'); rec_Neuron; source; visual; location; index];
NextObs = rec_Neuron.V';
IsDone = r < -30;

% IsDone is false, get reward
if ~IsDone
    Reward = 1;
else
    Reward = r;
end

end