% angle = getVisualAngle(visual,recon,source)
function angle = getVisualAngle(visual,recon, Nsources)
   if sum(visual.V_V)==0
      angle = 0;
   else
      if Nsources==1
         [~,locs] = max(visual.V_V);
         angle    = recon.neuron(locs);
      else
         [~,locs] = findpeaks(visual.V_V);
         if isempty(locs)
            angle = visual.V_V*recon.neuron; % / sum(jeff.V_A);
         else
            angle = recon.neuron(locs);
         end
      end
   end
end