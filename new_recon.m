function recon = new_recon(NR, type, gaussian_data, DEBUG, NOUT)
% The function recon_neuron = pl_recon(n_recon) returns the object
% recon_neuron, which represents the neurons responsible for the location
% detection of a sound. These neurons will be related to neurons generated
% from the pl_jeffress function and the visual input, that is, by calculating
% the difference between them. 
% n_recon is the number of reconstruction neurons, i.e. how many neurons
% are responsible for calculating the difference between visual and
% auditory sensor's potential.
% type represents the spacing of the neurons. The type 'u' represents
% uniform while 'g' represents gaussian. If gaussian is chosen, third 
% argument of the function input is needed, which represents gaussian data.
% Type in the mean as the first element and the variancec as the second
% element. 
% These neurons generated will always be from -pi/2 to pi/2.
   if type == 'l'
       recon.neuron = linspace(-pi/2, pi/2, NR);
       recon.neuron = recon.neuron'; % Transpose reconstruction  neuron
   elseif type == 'g'
       mu = gaussian_data(1);
       sigma = gaussian_data(2);
       start = normcdf(-pi/2,mu,sigma);
       finish = normcdf(pi/2,mu,sigma);
       increment = (finish-start)/(n_auditory-1);
       currentprob = start;
       for i = 1:NJ
           recon.neuron(i) = norminv(currentprob, mu, sigma);
           currentprob = currentprob + increment;
       end
       if recon.neuron(end)~=pi/2
           recon.neuron(end) = pi/2;
       end
       recon.neuron = recon.neuron';
   else
       display('type must be either ''u'' or ''g''');
       return
   end

   if DEBUG
      recon.DEBUG.err = zeros(NOUT, NR); 
   end
end