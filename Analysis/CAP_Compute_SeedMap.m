% tc: cell array with subject time courses
% seed: vector with seed information
function [SM,ASM] = CAP_Compute_SeedMap(tc,seed,seed_idx)

    tmp = logical(seed(:,seed_idx));

    % Subjectwise seed map, i.e. correlation between the average seed trace
    % and any other voxel. Size: n_voxels x 1
    if sum(tmp) == 1
        SM = cellfun(@(x) (corr(x,x(:,tmp))),tc,'un',0);
    else
        SM = cellfun(@(x) (corr(x,mean(x(:,tmp),2))),tc,'un',0);
    end

     % Average seed map across subjects from the population
    ASM = (mean(cell2mat(SM),2))';
end