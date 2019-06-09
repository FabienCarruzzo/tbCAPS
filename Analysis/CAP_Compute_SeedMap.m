% tc: cell array with subject time courses
% seed: vector with seed information
function [SM,ASM] = CAP_Compute_SeedMap(tc,seed,is_SS)

    % If the seed(s) entered are similar across subjects, then we just plot
    % the connectivity maps for the first one
    if ~is_SS
        tmp = logical(seed(:,1));

        % Subjectwise seed map, i.e. correlation between the average seed trace
        % and any other voxel. Size: n_voxels x 1
        if sum(tmp) == 1
            SM = cellfun(@(x) (corr(x,x(:,tmp))),tc,'un',0);
        else
            SM = cellfun(@(x) (corr(x,mean(x(:,tmp),2))),tc,'un',0);
        end

         % Average seed map across subjects from the population
        ASM = (mean(cell2mat(SM),2))';
        
    % Else, we compute FC with the proper seed every time
    else
        
        for s = 1:size(seed,2)
            seed_cell{s} = logical(seed(:,s));
        end
        
        SM = cellfun(@(x,y) (corr(x,mean(x(:,y),2))),tc,seed_cell,'un',0);
        
        ASM = (mean(cell2mat(SM),2))';
    end
end