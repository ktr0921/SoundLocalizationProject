function output = ddpg_Jeff(jeff, option)
% If struct2mat, then transform structure datatype to matrix
if strcmp(option,'struct2mat')
    % Transform jeff structure into 3th dim matrix (DEBUG is not stored)
    output = [];
    output(:, :, 1)   = jeff.J;
    output(:, 1, 2)   = jeff.neuron;
    output(:, 2, 2)   = [jeff.delta_lr; jeff.timeconst; 
        jeff.threshold; jeff.Llastspike; 
        jeff.Rlastspike; jeff.max_wgt; 
        jeff.min_wgt];
    output(:, 3, 2)   = jeff.V_L;
    output(:, 4, 2)   = jeff.V_R;
    output(:, 1:5, 3) = reshape(jeff.Lspikes, [20, 5]);
    output(:, 1:5, 4) = reshape(jeff.Rspikes, [20, 5]);
    output(:, :, 5)   = jeff.corr;
else % Else, transform matrix to structure datatype
    % Transform jeff 3d dim matrix into structure (DEBUG is not stored)
    output = struct();
    output.J = jeff(:, :, 1);
    output.neuron = jeff(:, 1, 2);
    output.delta_lr = jeff(1, 2, 2);
    output.timeconst = jeff(2, 2, 2);
    output.threshold = jeff(3, 2, 2);
    output.Llastspike = jeff(4, 2, 2);
    output.Rlastspike = jeff(5, 2, 2);
    output.max_wgt = jeff(6, 2, 2);
    output.min_wgt = jeff(7, 2, 2);
    output.V_L = jeff(:, 3, 2);
    output.V_R = jeff(:, 4, 2);
    output.Lspikes = reshape(jeff(:, 1:5, 3), [100, 1]);
    output.Rspikes = reshape(jeff(:, 1:5, 4), [100, 1]);
    output.corr   = jeff(:, :, 5);
end



