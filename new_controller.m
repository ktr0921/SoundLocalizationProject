% controller = new_controller(NJ, DEBUG)
%
function controller = new_controller(control_type, DEBUG, NOUT, varargin)
   NR = varargin{1};
   if strcmp('lqr', control_type)
      controller.xI = zeros(1, NR);
   elseif strcmp('lqr_ff', control_type)
      % Controller:   u = K(x - x_d) + u_d
      if nargin < 7, error('Controller lqr_ff needs 7 inputs'); end
      NA = varargin{2};
      eta = varargin{3};
      alpha = varargin{4};
      
      controller.x_d = rand(NA,NR);
      controller.u_d = rand(NA,NR);
      controller.u = rand(NA,NR);
%       controller.kn = new_feedback(NA, eta, alpha);
      controller.K = diag(rand(1,NA));
   end
   if DEBUG
      controller.DEBUG.u = zeros(NOUT, NR); 
   end
end