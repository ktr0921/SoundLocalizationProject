function mat_Out = ddpg_Transform(mat_In, option)

% [jeff, rec_Neuron, source, visual, location, index]
mat_Out = [];
if strcmp(option,'func2vec')
    % Storing jeff matrix in 4th dim matrix (DEBUG is not stored)
    mat_Out(:, :, 1, 1)   = mat_In(1).J;
    mat_Out(:, 1, 2, 1)   = mat_In(1).neuron;
    mat_Out(:, 2, 2, 1)   = [mat_In(1).delta_lr; mat_In(1).timeconst; 
        mat_In(1).threshold; mat_In(1).Llastspike; 
        mat_In(1).Rlastspike; mat_In(1).max_wgt; 
        mat_In(1).min_wgt];
    mat_Out(:, 3, 2, 1)   = mat_In(1).V_L;
    mat_Out(:, 4, 2, 1)   = mat_In(1).V_R;
    mat_Out(:, 1:5, 3, 1) = reshape(mat_In(1).Lspikes, [20, 5]);
    mat_Out(:, 1:5, 4, 1) = reshape(mat_In(1).Rspikes, [20, 5]);
    mat_Out(:, :, 5, 1)   = mat_In(1).corr;
    
    mat_Out(:, 1, 1, 2) = mat_In(2).neuron; % Store rec_Neuron
    mat_Out(:, 1, 1, 3) = mat_In(2).neuron; % Store rec_Neuron
    mat_Out(:, 1, 1, 4) = mat_In(2).neuron; % Store rec_Neuron
    mat_Out(:, 1, 1, 5) = mat_In(2).neuron; % Store rec_Neuron
    mat_Out(:, 1, 1, 6) = mat_In(2).neuron; % Store rec_Neuron
else
    mat_Out(:, :, 1) = mat_In(1).J;
    mat_Out(:, :, 2) = mat_In(1).J;
    mat_Out(:, :, 3) = mat_In(1).J;
    mat_Out(:, :, 4) = mat_In(1).J;
end

mat_In

jeff_J = jeff.J;
jeff_J_Dot = jeff.Jdot;
jeff_V_IA = jeff.V_IA;



