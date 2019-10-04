% I got angles by finding peaks, but for some reason this very occasionally 
% didn't work so well for only a single source, though mostly it was fine. 
% But seeing as I know it's a single source, just find the max.
T = true; F = false;
resetTime      = F;
resizeOutput   = F;
reformatOutput = F; 
% resizeWeights  = T;

plotAngles     = F; 
plotWeights    = T;

if resetTime
   NT = tstep;
   fin_time = NT*dt;
   t  = t(1:NT);
end
if resizeOutput
   % resize to ditch unsimulated time steps
   output_A_non = output_A_non(1:tstep,:);
   output_A_lqr = output_A_lqr(1:tstep,:);
   output_A_uns = output_A_uns(1:tstep,:);
   output_V     = output_V(1:tstep,:);
end

if reformatOutput
   % Franosch - non, unsupervised - uns, controlled - lqr
   plast  = {'vis','non','lqr'}; NP = length(plast);
   labels = {'Visual', 'Hebbian', 'Controlled', 'Unsupervised'};
   output = {output_V, output_A_non, output_A_lqr, output_A_uns};
   angles = cell(NP, 1);
   for ti=1:NT
      for pp=1:NP
         if ~any(output{pp}(ti,:))
            angles{pp}(ti) = NaN;
         else
            [~,loc]=(max(output{pp}(ti,:)));
            angles{pp}(ti) = recon_lqr.neuron(loc);
         end

      end
   end
end

if plotWeights
   f1 = figure; 
   imagesc(jeff_non.J'); %,[0 1]); axis image; colorbar;
   ylabel('Recon neurons');
   xlabel('Auditory neurons');
   set(f1,'name','Unsupervised Weights');
   
   f2 = figure; 
   imagesc(jeff_lqr.J'); %,[0 1]); axis image; colorbar;
   ylabel('Recon neurons');
   xlabel('Auditory neurons');
   set(f2,'name','Supervised Weights');
end

fopts = {'fontsize',24,'fontweight','bold'};

% Plotting the Angles (estimated vs real)
if plotAngles
   fh=figure; hold on; cols = ['k','b','r']; lw = [3 2 2];
   for pp=1:NP
      plot(t*1e3, angles{pp}/pi*180, cols(pp),'linewidth',lw(pp));
   end
   xlim([0 t(end)*1e3]);
   lh = legend(labels);
   set(gcf,'name','Angle estimates');
   xlabel('Time (s)',fopts{:});
   ylabel('Angle (degrees)',fopts{:});
   set(gca,fopts{:});
   set(lh,'box','off');
   saveFigure(fh,'AngleEstimate','fig','jpg');
end
