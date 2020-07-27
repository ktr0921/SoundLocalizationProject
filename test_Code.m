
clear;
clc;

configuration;
NOUT           = ceil(length(t)/out_dt);
[InitialObservation, LoggedSignals] = ddpg_Reset_Func();
control_non = new_controller('non', DEBUG, NOUT, NR);
Alpha = 1e-4*ones(NA,NR);

for i = 1:5000
    output_recon(i,:) = LoggedSignals.rec_Neuron.V;
    output_A_non(i,:) = LoggedSignals.jeff.V_A;   % auditory input into the recon neurons
    tstep = LoggedSignals.tstep;
    if mod(tstep, out_dt)==0
        OUT = tstep/out_dt;
    else
        OUT = false;
    end
    Action = controller(tstep, control_non, 'non', LoggedSignals.rec_Neuron, LoggedSignals.jeff, Alpha, dt, OUT);
    % disp(Action.u)
    [NextObs, Reward, IsDone, LoggedSignals] = ddpg_Step_Func(Action, LoggedSignals);
    
    % disp(Reward);
    % disp(IsDone);
end



%%
cdata = [45 60 32; 43 54 76; 32 94 68; 23 95 58];
xvalues = {'Small','Medium','Large'};
yvalues = {'Green','Red','Blue','Gray'};
h = heatmap(xvalues,yvalues,cdata);

h.Title = 'T-Shirt Orders';
h.XLabel = 'Sizes';
h.YLabel = 'Colors';
%%

test_List1 = test_List(1:100:end, :);

heatmap(test_List1)
h.Title = 'T-Shirt Orders';
h.XLabel = 'Sizes';
h.YLabel = 'Colors';
%%

X = 0;
DisplacementThreshold = 1;
Theta = 3;
AngleThreshold = 2;

IsDone = abs(X) > DisplacementThreshold || abs(Theta) > AngleThreshold;

% Get reward
if ~IsDone
    disp('Hi')
%     Reward = RewardForNotFalling;
else
    disp('No')
%     Reward = PenaltyForFalling;
end

disp(IsDone)

%%



NT     = length(t);
NOUT   = ceil(length(t)/out_dt);

jeff_non    = new_jeffress(NA, NR, 'l', 'stdp_fra', [], delta_lr, tau, neuronType, max_wgt, min_wgt, ...
    threshold, window_size, R_LIF, tau_m, v_reset, t_refrac, DEBUG, NOUT);
recon_non   = new_recon(NR, 'l', [], DEBUG, NOUT);
control_non = new_controller('non', DEBUG, NOUT, NR);
output_A_non= zeros(NT, NR);     % length(t) x NR
angle_A_non = cell(NT,1);        % auditory angle
% Create a visual (supervisor) neurons
visual         = new_supervisor(NV,'l',[],sigma_F, recon_non, neuronType, window_size, DEBUG, NOUT);

output_V       = zeros(NT, NR);     % length(t) x NR
angle_V        = cell(NT,1);        % visual (supervisor) angle
location       = cell(NT,1);        % source location

% initialise source data. There are Nsource sound sources, so
% initialise them one by one.
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

if mod(((tstep)/NT*100),10)==0
    fprintf('ITD simulation...%d%% complete...\n',uint16(round(tstep-1)/NT*100));
    save('non_lqr','*','-v7.3')
end

% if debugging & output time sample, set OUT to timestep
if mod(tstep, out_dt)==0
    OUT = tstep/out_dt;
else
    OUT = false;
end
tt = t(tstep);
% There are Nsource sources. Check each one for a particular time
% instance. Which one of the mastersource data is correct for that
% particular time instance?
for nsource = 1:Nsource
    % Extract the correct source data from the mastersource.
    % Mastersource variable is an array of struct. There are
    % Nsource elements of mastersource, where each element consists
    % of a struct.
    source(nsource)  = struct_source(tt, mastersource(nsource), index(nsource), source(nsource));
    % The infchange indicates that there is an information change
    % in the source. Therefore, this means that at this time, the
    % correct data that we need to extract from mastersource has
    % changed. Because it is already sorted such that the time is
    % increasing, then we can just simply increase the index for
    % that particular sound source.
    if source(nsource).infchange
        if doFranosch, jeff_non.corr = zeros(NA,NR);end  % Initialise the correlation term
        if doControl, jeff_ctl.corr = zeros(NA,NR);end
        if doInhibit, jeff_inh.corr = zeros(NA,NR);end
        index(nsource) = index(nsource) + 1;
    end
    location{tstep}(nsource) = source.location;
end
% Update the visual potential
visual      = update_supervisor(tstep, dt, visual, source, epsp, neuronType, OUT);

% Update the Jeffress ladder with the current auditory source
% Update the reconstruction neuron potential (calculate error in auditory estimate)
% Update the control signal
% Update auditory weights, J (from Jeffress ladder to recon neurons)
if doFranosch
    jeff_non    = update_jeffress(tstep, dt, jeff_non, source, epsp, neuronType, OUT, 1e-100);
    recon_non   = update_reconstruction(mode,recon_non,jeff_non,visual,inhibition,OUT);
end

%%
a = false;
if ~a
    disp(0)
end














