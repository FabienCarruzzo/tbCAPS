%% From the raw voxelwise data, computes a spatially averaged seed signal
% tcvox is a cell array (each cell has time x masked voxels dimension)
% seed is the seed mask (masked voxel x 1)
% S is a cell array (each cell is a time x 1 vector)
function [S] = CAP_get_seed_traces(tcvox,seed,SignMat)

    if SignMat(1)
        % Computation of the spatially averaged seed signal
        S = cellfun(@(x) zscore(mean(x(:,seed),2)), tcvox, 'un', 0);
    elseif SignMat(2)
        % Computation of the spatially averaged seed signal
        S = cellfun(@(x) (-1)*zscore(mean(x(:,seed),2)), tcvox, 'un', 0);
    else
        errordlg('PROBLEM WITH SIGN MAT !!!');
    end
end