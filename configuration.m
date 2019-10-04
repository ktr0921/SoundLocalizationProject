%% TO DO
% - multiple sources
% - jeffress ladder for realistic speed of sound, & multiple periods of the
%   frequency (i.e. don't assume that there will only be a single peak for
%   an input sound wave in the jeffress ladder, but periodic peaks)
% - try different controllers
% - parameter analysis, such as Q, R 

% Configuration parameters
clear source mastersource jeff visual controller recon

% General configt
DEBUG     = true;
init_time = 0;                      % intial time in second
fin_time  = 5e-3;                   % final time in second 
dt        = 1e-6;                   % sampling period in second per sample
ne_dist   = 1;                      % distance between cortical neuron
ne_vel    = 10;                     % velocity of signal travelling in cortical neuron
con_type  = 'lqr_ff';               % controller: 'non' (Franosch), 'lqr' (LQR), uns (unsupervised)
neuronType= 'r';                    % type 'r' for rate-based model and 's' for spike-based model
t         = (init_time:dt:(fin_time-dt))';
linear_loc= false;                  % make location linearly increase for testing
out_dt    = 1e0;                    % time resolution of debug output
inhibition = false;                 % lateral inhibition of reconstruction neurons

% Jeffress, supervisor and reconstruction neurons config
NA = 20;                         % # of auditory nuerons
NV = 20;                         % # of visual neurons
NR = 20;                         % # of reconstruction neurons

% Source config
Ntrials   = 100;
Nsource   = 1;
t_src     = sort([zeros(Nsource, 1) rand(Nsource,Ntrials-1)],2) * fin_time; % sort ascending along time dimension
freq = ones(size(t_src)) * 5e3;
tol       = 0.4;         % delta_lr has to be strictly less than 1/freq; downscale to tol/freq to ensure this condition
delta_lr  = tol/freq(1);  % the time taken for the signal to reach right end from the left end of the Jeffress ladder and vice versa
Vmax      = 10;
ampl_A    = ones(size(t_src)) * Vmax / (NA);
ampl_V    = ones(size(t_src)) * Vmax;
% lctn      = rand(size(t_src))*pi - pi/2;  % randomize between -pi/2 and pi/2
lctn      = linspace(-pi/2,pi/2,length(t_src)); % Hardcode angles (comment this line and uncomment the previous one for random angle)
lctn_std = 0;                    % Gaussian std dev for VISUAL input activity in spikes / sec (How the potential is distributed along the visual neurons)
base     = 1e-1;                    % ampl_A;  % baseline in spikes / sec (to prevent negative rate)
sigma    = 0;                    % Gaussian standard deviation for white noise in spikes / sec

sigma_F  = 0.25;                 % std dev of visual weights, F_ij (if 0 -> F is identity matrix)
tau      = 1e1;                 % learning time constant, tau (reciprocal of eta in Franosch)
epsp_str = 'expon';              % 'delta' for Dirac Delta or 'expon' for hyperbolic exponential
t_r      = 1e-6;                 % rise time for exp fn (for epsp)
max_wgt  = 1;
min_wgt  = 0;

% neuron specific configuration parameters
switch neuronType
   case 's'
        % Parameters for Leaky Integrate-and-Fire (LIF)
        R_LIF     = 1e+3;               % The 'resistance' of LIF (leaky term)
        tau_m     = 1e-6;               % The time constant of the LIF equation
        t_refrac  = 1e-4;               % refractory period for LIF neurons
        v_reset   = 0;                  % The value of potential the neuron will reset to after spike
        threshold = 0.04; % R_LIF*(ampl_A(1)+base(1))*t_r*0.99;

        % window_size is used for convolution later - need in timsteps not time.
        window_width  = t_refrac; % need to recall spikes & potential for refrac period
        window_size = round(min([window_width/dt fin_time/dt]));
   case 'r'
        threshold = 0.9;
%         threshold = 0.1;
        R_LIF = [];
        tau_m = [];
        t_refrac  = 1e-4; % refractory period for LIF neurons
        v_reset = [];
        
        % Initialize variables that are useful only for 's'  
        window_width  = t_refrac; % need to recall spikes & potential for refrac period
        window_size = round(min([window_width/dt fin_time/dt]));
end


if DEBUG && linear_loc
   % Set min time for each angle to be 3 cycles of a sound wave
   t_min = (1./freq)*3;
   if sum(t_min)>fin_time
      fin_time = sum(t_min);
      warn     = sprintf('ITD:inputWarning: less than 3 cycles per sound source\n\n');
      cprintf('*Blue',sprintf('\nITD:inputWarning: less than 3 cycles per sound source\n\n'));
      return;
   end
%    lctn  = linspace(-pi/2, pi/2, Ntrials);
   lctn = linspace(-pi/2,pi/2,length(t_src)); % Hardcode angles
end

% The area under the epsp function must be 1. Therefore, in discrete form,
% this has to be equal to the sampling period. Refer to the book for the
% reason why.
if epsp_str == 'delta'
    epsp = 1.0/dt*[1 zeros(1,length(t)-1)];
elseif epsp_str == 'expon'
    epsp = t(2:end)/t_r.*exp(-t(2:end)/t_r);
    epsp = epsp/(dt*sum(epsp));
end

if strcmp('s',neuronType)
    epsp = epsp(1:min([window_width/dt fin_time/dt])); % This line works with spike based
end

mode  = 'sup';       % the mode of the system, either unsupervised ('uns') or supervised ('sup')
model = 'stdp_var';  % the model of the system ('stdp_var' for variation of stdp and 'stdp_fra' for Franosch stdp formula)
Alpha = 1e-4*ones(NA,NR); % the rate of change of J depends on itself for stdp_var model. Choose this to be constant.
% NOTE: Alpha represents the variable gamma in Franosch stdp. Alpha(i,j) is
% gamma responsible for J(i,j).

% Construct the struct of the source config to a variable mastersource
for nsource = 1:Nsource
    mastersource(nsource).t_source  = t_src(nsource,:);
    mastersource(nsource).frequency = freq(nsource,:);
    mastersource(nsource).amplitude = ampl_A(nsource,:);
    mastersource(nsource).visualAmp = ampl_V(nsource,:);
    mastersource(nsource).location  = lctn(nsource,:);
    mastersource(nsource).loc_std   = lctn_std;
    mastersource(nsource).baseline  = base;
    mastersource(nsource).gauss_std = sigma;
end

% PARSING ALL THE INPUTS
% Parsing init_time
if ~isnumeric(init_time)
    error('init_time has to be a number');
    return
elseif init_time<0
    error('init_time has to be non-negative number');
    return
end

% Parsing fin_time

% Parsing samp_per

% Parsing window_width

% Parsing ne_dist

% Parsing ne_vel

% Parsing con_type