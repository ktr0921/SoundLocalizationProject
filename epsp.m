% [potential, kernel] = epsp(v,T,dt,kernelType,params,doplot)
%
% Generate an excitatory post synaptic potential kernel. Convolve with
% neuron potential, v, if v is not empty or simply equal to 1. Else return
% the kernel window only.
% Inputs:
%   v           - neuron voltgage potential to convolve with the kernel
%                 (set to 1 to ignore & return the kernel window instead)
%   T           - time length of the kernel
%   dt          - time resolution of the kernel
%
%   kernelType  - specify epsp convolution kernel
%       'gaussian': a gaussian kernel, specified by sigma (default)
%       'poisson':  poisson kernel 
%       'exp':      a combination of exponential growth and decay funcions 
%                   that resemble a postsynaptic potential.
%   params  - input parameters required for the chosen kernel type
%       gaussian:   sigma (defaults to 5 samples)
%       poisson:    no input parameters
%       exp:        [tau_d tau_r] to specify decay & rise (opt.) time 
%                   cosntants (defaults to 10 and 0, wrt Thompson et al. 1996)
%
% Outputs:
%   potential   - neuron output potential, derived from spikes convolved
%                 with the kernel window
%   window      - window generated from the kernel type requested
function [potential, window] = epsp(varargin)
  %% Parse inputs
    v_      = 1;
    dt_     = 0.001;
    T_      = 10;
    kernel_ = 'exp';
    params_ = [];
    doplot_ = false;
    
    Tr_     = []; % 1e-2;
    Td_     = 1e-2;
    sigma_  = 1e-2;
    lambda_ = 5;
    
    nargs   = length(varargin);
    optargs = {v_, [], [], kernel_, params_, doplot_}; % T, dt, kernel, params, doplot
    optargs(1:nargs) = varargin;
    [v, T, dt, kernel, params, doplot] = optargs{:};
    
    if isempty(v),      v      = v_;        end
    % If voltage is a vector, then we only need one of T, dt to determine
    % the time vector
    if isempty(T) && (~isempty(dt) && ~isscalar(v))
        T = length(v)*dt;
    elseif isempty(dt) && (~isempty(T) && ~isscalar(v))
        dt = T/length(v);
    end
    if isempty(T),      T      = T_;        end
    if isempty(dt),     dt     = dt_;       end
    if isempty(kernel), kernel = kernel_;   end
    if isempty(doplot), doplot = doplot_;   end
    if ~isfloat(v),     v      = double(v); end
    [nT, nN] = size(v); % number of time points, & number of neurons
    
    % If unknown kernel type, set to default kernel (gotta do it up here so
    % that default parameter inputs are correctly set)
    switch lower(kernel)
        case {'gaussian','gauss','g','epsp','eps','psp','e','poisson','p','poiss','pois'}
        otherwise
            kernel = kernel_;
    end
    % Set default kernel parameters if none supplied
    if isempty(params)
        switch lower(kernel)
        case {'gaussian','gauss','gaus','g'}
                params = sigma_;
        case {'exp','decap','exponential','epsp','epse','e','esp','eps','psp'}
                params = [Tr_, Td_];
            case {'poisson','p','poiss','pois'}
                params = [lambda_];
        end
    end
    % Populate relevant kernel parameters from param input
    switch lower(kernel)
        case {'gaussian','gauss','gaus','g'}
            sigma = params(1);
        case {'exp','decay','exponential','epsp','epse','e','esp','eps','psp'}
            if length(params)>=2
                Td = params(1);
                Tr = params(2);
            elseif isscalar(params) % only input 1 of the exp values (Te)
                Td = params;
                Tr = Tr_;
            else
                Td = Td_;
                Tr = Tr_;
            end
        case {'poisson','p','poiss','pois'}
            lambda = params;
        otherwise
            kernel = kernel_;
    end
    
  %% Build kernel
    time     = (dt:dt:T)';
    nsamples = length(time);

    switch lower(kernel)
        case {'gaussian','gauss','gaus','g'}
            window  = normpdf(time-mean(time),0,sigma);
        case {'poisson','poiss','pois','p'}
            window  = poisspdf(time/dt,lambda);
        case {'exp','decap','exponential','epsp','epse','e','esp','eps','psp'}
            if isempty(Tr) || Tr==0
                window = exp(-time/Td); % normalised so coeff 1/Td irrelevant
            else
                window = (1-exp(-time/Tr)).*(exp(-time/Td));
                window  = window./sum(window);
            end
    end
    
%% Convolve
    % if input voltage is just 1, only the kernel is required. O/wise, we
    % need to convolve the epsp kernel with the neuron's voltage potential
    if isempty(v) || isscalar(v)
        calcPotent = false;
%         potential  = [];
        potential  = window;
        vtime      = time;
    else
        calcPotent = true;
        vtime      = (1:length(v))*dt;
        if isDim(v,2)
            % convolve with 'same' shifts everything in the wrong direction
            % (makes the convolved action earlier by length(window) time
            % points). Faaaark! Instead, use 'full' & cull the last
            % length(window)-1 time points to make the size match v
            potential = arrayfun(@(vi) conv(v(:,vi),window,'full'), ...
                                 1:nN, 'uniformoutput', false);
            potential = cell2mat(potential);
            potential = potential(1:end-length(window)+1,:);
        else
            potential = conv(v,window,'same'); % output same length as input voltage
        end
    end

    if doplot
        if calcPotent
            nr = 2; nc = 1; si=1;
        else
            nr = 1; nc = 1; si=1;
        end
        figure;
        subplot(nr,nc,si), si=si+1;
        plot(time,window);
        title(sprintf('%s epsp kernel',upper(kernel)));

        if calcPotent
            subplot(nr,nc,si), si=si+1;
            plot(vtime,potential)
            title('Neuron potential');
        end
    end
end





