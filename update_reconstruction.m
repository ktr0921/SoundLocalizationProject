function recon = update_reconstruction(mode,recon,jeff,visual,inhibition,DEBUG)
% the function recon = update_reconstruction(recon,jeffress,supervisor)
% updates the potential inside the reconstruction neurons recon. The input
% is the output potential of Jeffress Ladder 'jeffress' and the output
% potential of the supervisor 'supervisor'. The potential from Jeffress
% will be weighted by J_ij correspondingly to create the final potential
% from the slave neuron, while the potential from the visual system will be
% weighted by F_pq correspondingly to create the final potential from the
% supervisor neurons. These final potentials will then be compared to get
% the updated potential for the reconstruction neuron, W.
   if mode=='sup'
      V_visual = visual.V_V;
      V_audio  = jeff.V_A;
      err      = - V_visual + V_audio;
      recon.V  = err;% .* heaviside(err);
      if inhibition
        recon.V = get_inhibition(recon.V);
      end
   elseif mode=='uns'
      V_visual = visual.V_V;
      V_audio  = jeff.V_A;
      recon.V  = V_visual + V_audio;
      err = 0;
   end
   if DEBUG
      recon.DEBUG.err(DEBUG,:) = err;
   end
end

function v = get_inhibition(v)
    % Inhibition of neighbours
    inhibit = v;
    inhibit(inhibit < 1)= 1;
    inhibit_ = ones(size(inhibit));
    for i = 1:length(v)
      if i > 1 && i < length(v)
        inhibit_(i) = 2./(inhibit(i - 1) + inhibit(i + 1));
      elseif i == 1
        inhibit_(i) = 1./(inhibit(i + 1));
      else
        inhibit_(i) = 1./(inhibit(i - 1));
      end
    end
    v = v .* inhibit_;
end