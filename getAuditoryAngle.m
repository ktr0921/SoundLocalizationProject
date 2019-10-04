% angle = getAuditoryAngle(jeff,recon)
function angle = getAuditoryAngle(jeff,recon, Nsources)
   if length(unique(jeff.V_A))==1
      angle = NaN;
   elseif Nsources==1
      [~,locs] = max(jeff.V_A);
      angle = recon.neuron(locs);
   else
      [~,locs] = findpeaks(jeff.V_A);
      if isempty(locs)
         angle = jeff.V_A*recon.neuron; % / sum(jeff.V_A);
      else
         angle = recon.neuron(locs);
      end
   end
end