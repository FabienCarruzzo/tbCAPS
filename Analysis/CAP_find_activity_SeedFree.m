%% Finds the moments of (de)activation in a group of fMRI subjects
% Inputs
% tcvox: cell aray with a seed signal in each cell (time points x masked voxels)
% seed: masks for the seeds used (masked voxels x n_seed)
% T: threshold for the retention of active or inactive frames
% FDall: framewise displacement traces for all subjects (time points x n_subjects)
% Mot_thresh: threshold (in mm) to use for scrubbing
%
% Outputs
% Xonp: cell array (each cell dimension masked voxels x n_retained_frames)
% for active frames brain patterns
% p: 5xn_subject matrix with percentages of scrubbed frames, scrubbed and
% active frames, scrubbed and inactive frames, retained frames for active
% analysis and retained frames for inactive analysis
function [Xonp,p,Indices,Xonp_scrub] = CAP_find_activity_SeedFree(tcvox,FDall,Mot_thresh)

    % Computes the indices of the peak values that must be excluded
    flag = FDall>Mot_thresh;
    
    % We want to store the indices of scrubbed frames and return it later
    Indices.scrubbed = logical(flag);
    
    % Each cell of flag will contain a time x 1 vector of logicals (1 if
    % the frame is to be censored, 0 otherwise)
    flag = num2cell(flag,1);
    
    % 1 x n_subject vector containing the percentage of scrubbed frames (throughout the whole scan) 
    p_scrubbed = cell2mat(cellfun(@(x) sum(x)/length(x)*100, flag,'un',0));
    
    % We will want to retain all frames (prior motion scrubbing)
    xindp = {};
    for s = 1:size(FDall,2)
        xindp{s} = logical(ones(size(FDall,1),1));
    end
     
    % flag now contains the traces with high activity AND high motion
    flag_active = cellfun(@(x,y) x & y,xindp, flag,'un',0);

    % Vector (1xn_subj) with the percentage of traces removed because of too high
    % motion and being selected as active
    p_scrubactive = cell2mat(cellfun(@(x) sum(x)/length(x)*100, flag_active,'un',0));

    % My indices of active/inactive frames now contain only the non
    % corrupted frames
    xindp = cellfun(@(x,y) x & ~y, xindp,flag_active,'un',0);

    Indices.kept.active = cell2mat(xindp);
    Indices.scrubbedandactive = cell2mat(flag_active);
    
    % Each cell contains the frames selected as active or as inactive (if
    % inactive, the sign is reversed, i.e. inactivation is a positive
    % signal). Size: masked voxels x n_retained_frames
    Xonp = cellfun(@(x,y) x(y,:)',tcvox,num2cell(Indices.kept.active,1),'un',0);
    Xonp_scrub = cellfun(@(x,y) x(y,:)',tcvox,num2cell(logical(Indices.scrubbedandactive),1),'un',0);
    
    % Percentage of active and inactive frames retained per subject
    p_active = cell2mat(cellfun(@(x) size(x,2)/size(FDall,1)*100, Xonp,'un',0));
    
    % Matrix containing all the interesting probabilities
    p = [p_scrubbed; p_scrubactive; p_active];
end