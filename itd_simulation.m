%% Questions
% - Change V_IA to be sin * epsilon, rather than just sin, to be the same
% as Franosch
% - Jdot - Ax is very small (-0.0167) while jeff.V_IA*control(loop) is on the
%     order of hundreds, meaning Ax is pretty irrelevant??
% - The weights for each auditory neuron seem to be doing the same thing?
% - Can you prove convergence (with a time constant?!) and stability?
% - Do you still have a work computer - you can run simulations on it?

%% TO DO:
% - find peaks in angle estimates that allow for interpolation btwn angles
% - When NR > NJ, every second weight is 0
% - multiple sources
% - jeffress ladder for realistic speed of sound, & multiple periods of the
%   frequency (i.e. don't assume that there will only be a single peak for
%   an input sound wave in the jeffress ladder, but periodic peaks)
% - try different controllers
% - parameter analysis, such as Q, R

% close all;
clear
% Run Simulation: True or False
T = true; F = false;
continueSim  = F;   % allows you to continue a simulation from where you left off
runSim       = T;
plotOutput   = T;   % Plot estimated auditory location vector & visual output
plotAngles   = T;   % Plot estimated and reference angles
plotWeights  = T;   % Plot Jeffress weight: True or False
plotJeffLIF  = F;   % Plot jeffress LIF neuron potential
plotAudInput = T;   % Plot Jeffress input (before weighting)

doFranosch   = T;
doControl    = T;
doUnsup      = F;
doInhibit    = inhibition;   % Lateral inhibition - defined in configuration

if runSim
    if continueSim
        startTime = tstep;
    else
        clear jeff recon visual source control
        configuration;
        tstart = tic;
        NT     = length(t);
        NOUT   = ceil(length(t)/out_dt);
        out_t  = (0:out_dt:(NT-1))'*dt;
        
        % Create a new Jeffress Ladder
        % Create reconstruction neurons
        % New controller
        % Initialise all variables which we want to monitor
        if doFranosch
            jeff_non    = new_jeffress(NA, NR, 'l', 'stdp_fra', [], delta_lr, tau, neuronType, max_wgt, min_wgt, ...
                threshold, window_size, R_LIF, tau_m, v_reset, t_refrac, DEBUG, NOUT);
            recon_non   = new_recon(NR, 'l', [], DEBUG, NOUT);
            control_non = new_controller('non', DEBUG, NOUT, NR);
            output_A_non= zeros(NT, NR);     % length(t) x NR
            angle_A_non = cell(NT,1);        % auditory angle
           % Create a visual (supervisor) neurons
           visual         = new_supervisor(NV,'l',[],sigma_F, recon_non, neuronType, window_size, DEBUG, NOUT);
        end
        
        if doControl
            jeff_ctl    = new_jeffress(NA, NR, 'l', 'stdp_fra', [], delta_lr, tau, neuronType, max_wgt, min_wgt, ...
                0.1, window_size, R_LIF, tau_m, v_reset, t_refrac, DEBUG, NOUT);
            recon_ctl   = new_recon(NR, 'l', [], DEBUG, NOUT);
            control_ctl = new_controller('lqr_ff', DEBUG, NOUT,NR,NA, 1/tau, Alpha);
            output_A_ctl= zeros(NT, NR);     % length(t) x NR
            angle_A_ctl = cell(NT,1);        % auditory angle
           % Create a visual (supervisor) neurons
           visual         = new_supervisor(NV,'l',[],sigma_F, recon_ctl, neuronType, window_size, DEBUG, NOUT);
        end
        
        if doInhibit
            jeff_inh    = new_jeffress(NA, NR, 'l', 'stdp_fra', [], delta_lr, tau, neuronType, max_wgt, min_wgt, ...
                0.1, window_size, R_LIF, tau_m, v_reset, t_refrac, DEBUG, NOUT);
            recon_inh   = new_recon(NR, 'l', [], DEBUG, NOUT);
            control_inh = new_controller('lqr_ff', DEBUG, NOUT,NR,NA, 1/tau, Alpha);
            output_A_inh= zeros(NT, NR);     % length(t) x NR
            angle_A_inh = cell(NT,1);        % auditory angle
           % Create a visual (supervisor) neurons
           visual         = new_supervisor(NV,'l',[],sigma_F, recon_inh, neuronType, window_size, DEBUG, NOUT);
        end
        
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
        
        startTime = 1;
        disp(['Running ''' neuronType ''' based']);
    end
    for tstep = startTime:NT-1
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
            control_non = controller(tstep, control_non, 'non', recon_non, jeff_non, Alpha, dt, OUT);
            jeff_non    = update_weight(tstep,jeff_non,'sup','stdp_var',source, control_non, dt, Alpha, OUT);
            output_recon(tstep,:) = recon_non.V;
            output_A_non(tstep,:) = jeff_non.V_A;   % auditory input into the recon neurons
            angle_A_non{tstep}    = getAuditoryAngle(jeff_non,recon_non,Nsource); % auditory angle est
            output_V(tstep,:)= visual.V_V; % visual input into the recon neurons
            angle_V{tstep}   = getVisualAngle(visual,recon_non,Nsource); % visual angle
        end
        
        if doControl
            jeff_ctl    = update_jeffress(tstep, dt, jeff_ctl, source, epsp, neuronType, OUT, 0.01);
            recon_ctl   = update_reconstruction(mode,recon_ctl,jeff_ctl,visual,inhibition,OUT);
            control_ctl = controller(tstep, control_ctl, 'lqr_ff', recon_ctl, visual ,jeff_ctl, Alpha, dt, OUT);
            jeff_ctl    = update_weight(tstep,jeff_ctl,'sup','stdp_var',source, control_ctl, dt, Alpha, OUT);
            output_recon(tstep,:) = recon_ctl.V;
            output_A_ctl(tstep,:) = jeff_ctl.V_A;   % auditory input into the recon neurons
            angle_A_ctl{tstep}    = getAuditoryAngle(jeff_ctl,recon_ctl,Nsource); % auditory angle est
            output_V(tstep,:)= visual.V_V; % visual input into the recon neurons
            angle_V{tstep}   = getVisualAngle(visual,recon_ctl,Nsource); % visual angle
        end
        
        if doInhibit
            jeff_inh    = update_jeffress(tstep, dt, jeff_inh, source, epsp, neuronType, OUT, 0.01);
            recon_inh   = update_reconstruction(mode,recon_inh,jeff_inh,visual,true,OUT);
            control_inh = controller(tstep, control_inh, 'lqr_ff', recon_inh, visual ,jeff_inh, Alpha, dt, OUT);
            jeff_inh    = update_weight(tstep,jeff_inh,'sup','stdp_var',source, control_inh, dt, Alpha, OUT);
            output_recon(tstep,:) = recon_inh.V;
            output_A_inh(tstep,:) = jeff_inh.V_A;   % auditory input into the recon neurons
            angle_A_inh{tstep}    = getAuditoryAngle(jeff_inh,recon_inh,Nsource); % auditory angle est
            output_V(tstep,:)= visual.V_V; % visual input into the recon neurons
            angle_V{tstep}   = getVisualAngle(visual,recon_inh,Nsource); % visual angle
        end        
        
        if mod( tstep, 500 ) == 0
            nan;
        end
    end
    tstop = toc(tstart);
    printTime(tstop,'Learning sound localisation took');
    save(sprintf('ITD_%s_%s.mat'));
end

%%
% Plotting the membrane potential (auditory vs visual)
if plotOutput
    if doFranosch, ff = figure; end
    if doControl, fc = figure; end
    fv = figure;
    fr = figure;
    nN   = min([30 NR]); % check first, last, & an even spread in between
    nr   = ceil(sqrt(nN)); nc = ceil(nN/nr);
    dN   = ceil(min([NA/(nN-1) 1]));
    nind = round(1:dN:NA); nind(end+1) = NA; nind(nind>NA)=[]; nind=unique(nind);
    
    %%
    for ii=1:nN
        if doFranosch
            figure(ff);
            subplot(nr,nc,ii), % subplot(nr,nc,switchRowsCols(nr,nc,ii)),
            plot(t, output_A_non(:,nind(ii)));
            ylabel(sprintf('%d [%.2g^o]',nind(ii), ((nind(ii)-1)/(NR-1)*pi - pi/2)/pi*180));
            xlim([init_time, fin_time]);
        end
        if doControl
            figure(fc);
            subplot(nr,nc,ii), % subplot(nr,nc,switchRowsCols(nr,nc,ii)),
            plot(t, output_A_ctl(:,nind(ii)));
            ylabel(sprintf('%d [%.2g^o]',nind(ii), ((nind(ii)-1)/(NR-1)*pi - pi/2)/pi*180));
            xlim([init_time, fin_time]);
        end
        figure(fv);
        subplot(nr,nc,ii), % subplot(nr,nc,switchRowsCols(nr,nc,ii)),
        plot(t, output_V(:,nind(ii)));
        ylabel(sprintf('%d [%.2g^o]',nind(ii), ((nind(ii)-1)/(NR-1)*pi - pi/2)/pi*180));
        xlim([init_time, fin_time]);
        % Plot reconstruction
        figure(fr);
        subplot(nr,nc,ii), % subplot(nr,nc,switchRowsCols(nr,nc,ii)),
        plot(t(2:end), output_recon(:,nind(ii)));
        ylabel(sprintf('%d [%.2g^o]',nind(ii), ((nind(ii)-1)/(NR-1)*pi - pi/2)/pi*180));
        xlim([init_time, fin_time]);
    end
    if doFranosch, set(ff,'name','Franosch auditory output'); end
    if doControl,  set(fc,'name','Controlled auditory output'); end
    if doUnsup,    set(fu,'name','Unsupervised auditory output'); end
    set(fv,'name','Visual output');
    set(fr,'name','Reconstruction output');
end

%% Plotting the Angles (estimated vs real)
if plotAngles
    figure; hold on; leg = []; plotLines = [];
    for ti=1:tstep
        lra = plot(ti*dt*ones(size(angle_V{ti})),angle_V{ti},'bo'); hold on;
        if ti==1, plotLines(1) = lra; leg{1} = 'Visual angle'; end
        if doFranosch
            lyf = plot(ti*dt*ones(size(angle_A_non{ti})),angle_A_non{ti},'ro');
            if ti==1, plotLines(end+1) = lyf; leg{end+1} = 'Franosch audio angle'; end
        end
        if doControl
            lyc = plot(ti*dt*ones(size(angle_A_ctl{ti})),angle_A_ctl{ti},'ko');
            if ti==1, plotLines(end+1) = lyc; leg{end+1} = 'Controlled audio angle'; end
        end
        if doInhibit
            lyi = plot(ti*dt*ones(size(angle_A_inh{ti})),angle_A_inh{ti},'gx', 'MarkerSize', 12);
            if ti==1, plotLines(end+1) = lyi; leg{end+1} = 'Inhibition audio angle'; end
        end
    end
    legend(plotLines, leg);
    set(gcf,'name','Angle estimates');
end
%%
if plotWeights % format of weights, J, is [time, NA, NR)
    num_plots = min([NA 20]);
    delta_n = floor(NA/num_plots); % round down so indices never go over NA
    n_ind   = 1:delta_n:NA;        % get neuron indices from throughout NA
    nr = ceil(sqrt(num_plots));
    nc = ceil(num_plots/nr);
    if doFranosch
        f2f = figure;
        imagesc(jeff_non.J',[0 1]); axis image; colorbar;
        ylabel('Recon neurons');
        xlabel('Auditory neurons');
        set(f2f,'name','Franosch weights');
    end
    if doControl
        f2c = figure;
        imagesc(jeff_ctl.J',[0 1]); axis image; colorbar;
        ylabel('Recon neurons');
        xlabel('Auditory neurons');
        set(f2c,'name','Control weights');
    end
    if doInhibit
        f2i = figure;
        imagesc(jeff_inh.J',[0 1]); axis image; colorbar;
        ylabel('Recon neurons');
        xlabel('Auditory neurons');
        set(f2i,'name','Inhibition weights');
    end
end

if plotAudInput
    if doFranosch, ff = figure; end
    if doControl,  fc = figure; end
    if doUnsup,    fu = figure; end
    nN   = min([30 NA]); % check first, last, & an even spread in between
    nr   = ceil(sqrt(nN)*2); nc = ceil(nN/nr);
    dN   = ceil(min([NA/(nN-1) 1]));
    nind = round(1:dN:NA); nind(end+1) = NA; nind(nind>NA)=[]; nind=unique(nind);
    if doFranosch, subplot(ceil(sqrt(nN)),ceil(sqrt(nN)),1), plot(t_src, lctn/pi*180, 'ro'); end
    xlim([init_time, fin_time]); ylabel('Angle'); % xlabel('Time');
    for ii=1:nN
        if doFranosch
            figure(ff);
            subplot(ceil(sqrt(nN)),ceil(sqrt(nN)),ii+1), % subplot(nr,nc,switchRowsCols(nr,nc,ii+1)),
            plot(out_t,jeff_non.DEBUG.V_IA(:,nind(ii)))
            ylabel(sprintf('%d [%.2g^o]',nind(ii), ((nind(ii)-1)/(NA-1)*pi - pi/2)/pi*180));
            %           ylim([0 2*max(ampl_A)]);
%             ylim([0 1]);
        end
        if doControl
            figure(fc);
            subplot(ceil(sqrt(nN)),ceil(sqrt(nN)),ii+1), % subplot(nr,nc,switchRowsCols(nr,nc,ii+1)),
            plot(out_t,jeff_ctl.DEBUG.V_IA(:,nind(ii)))
            ylabel(sprintf('%d [%.2g^o]',nind(ii), ((nind(ii)-1)/(NA-1)*pi - pi/2)/pi*180));
            %           ylim([0 2*max(ampl_A)]);
%             ylim([0 1]);
        end
    end
    if doFranosch, set(ff,'name','V_IA');end
    if doControl, set(fc,'name','V_IA');end
end

if plotJeffLIF
    flif = figure;
    nN   = 30; % check first, last, & an even spread in between
    nr   = ceil(sqrt(nN)*2); nc = ceil(nN/nr);
    dN   = min([NA/(nN-1) 1]);
    nind = round(1:dN:NA); nind(end+1) = NA; nind(nind>NA)=[];
    subplot(nr,nc,1), plot(t_src, lctn/pi*180, 'ro');
    xlim([init_time, fin_time]); ylabel('Angle'); % xlabel('Time');
    for ii=1:nN
        subplot(nr,nc,ii+1), % subplot(nr,nc,ii+1), % subplot(nr,nc,switchRowsCols(nr,nc,ii+1)),
        plot(out_t,jeff.DEBUG.V_LIF(:,nind(ii)))
        ylabel(sprintf('%d [%.2g^o]',nind(ii), ((nind(ii)-1)/(NA-1)*pi - pi/2)/pi*180));
        ylim([0 threshold]);
    end
    set(flif,'name','V_LIF');
end



