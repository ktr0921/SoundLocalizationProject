function source = struct_source(tt, mastersource, index, source)
% Construct a struct variable (source) given the current time tt, the
% mastersource main data and the index that the previous source sample is
% based on (see the index of mastersource). This function always looks ahead
% on the next mastersource.t_source. When the current time tt is greater
% than or equal to the next mastersource.t_source, it will change the current sample
% source information into the information provided in mastersource of the 
% corresponding next index. This function also returns the variable 
% infchange, telling us whether the sample source information is changed 
% after this function is called.

   % if there are no inputs, return an empty structure (for preallocation
   % of memory)
    if nargin==0
      source = struct('frequency',[],'amplitude',[],'visualAmp',[],'location',[],...
                      'loc_std',[],'baseline',[],'gauss_std',[],'t_source',[],'infchange',[]);
   	return;
    end

   t_src = mastersource.t_source; % time vector of sound changes
   max_index = length(t_src);
   if index <= max_index
      if index==1 % the sound has not begun yet
        source.infchange = false;      % no change induced in sound
      end
      % If current time is greater than the finish time of this sound, 
      % change it if there's another sound to change to
      if length(t_src)>=(index+1) && tt > t_src(index+1) 
        source.infchange = true;       % true if the information has changed        
      else
        source.infchange = false;      % no change in the sound yet
      end
      source.frequency = getVal(mastersource.frequency,index); 
      source.amplitude = getVal(mastersource.amplitude,index);
      source.visualAmp = getVal(mastersource.visualAmp,index);
      source.location  = mastersource.location(index);
      source.loc_std   = getVal(mastersource.loc_std,index);
      source.baseline  = getVal(mastersource.baseline,index);
      source.gauss_std = getVal(mastersource.gauss_std,index);
      source.t_source  = getVal(mastersource.t_source,index);
   end
end

% Allows user to input a vector of values, one for each trial, or a single
% value if they want it to be the same for all trials. 
function y = getVal(x,index)
   if isscalar(x), y = x; return; end
   y = x(index); return; 
end