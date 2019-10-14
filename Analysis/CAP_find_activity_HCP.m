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
function [Xonp,p,Indices,idx_sep_seeds,flag_outofcuriosity,Xonp_full] = CAP_find_activity_HCP(tcvox,seed,T,FDall,Mot_thresh,SelMode,SeedType,SignMat,is_SS)

    % Contains the indices (logical) of the selected frames for each seed
    % subcase
    idx_sep_seeds = nan(size(FDall,1),size(FDall,2),size(seed,2));

    % Computes the indices of the peak values that must be excluded. If we
    % go above threshold for motion, we censor the current frame, but also
    % the one just before and the one just after
    flag = FDall>Mot_thresh;
    flag2 = circshift(flag,[1 0]);
    flag3 = circshift(flag,[-1 0]);
    flag4 = circshift(flag,[-2 0]);
    flag5 = circshift(flag,[-3 0]);
    flag6 = circshift(flag,[-4 0]);
    flag7 = circshift(flag,[2 0]);
    flag8 = circshift(flag,[3 0]);
    flag9 = circshift(flag,[4 0]);
    flag10 = circshift(flag,[5 0]);
    flag11 = circshift(flag,[6 0]);
    flag12 = circshift(flag,[7 0]);
    flag = flag+flag2+flag3+flag4+flag5+flag6+flag7+flag8+flag9+flag10+flag11+flag12; 
    flag_outofcuriosity = flag;
    flag(flag>1)=1;
    
    % We want to store the indices of scrubbed frames and return it later
    Indices.scrubbed = logical(flag);
    
    % Same for the retained frames
    switch SeedType
        case 'Intersection'
            Indices.kept.active = logical(zeros(size(FDall,1),size(FDall,2)));
        otherwise
            Indices.kept.active = logical(ones(size(FDall,1),size(FDall,2)));
    end
    
    
    % Each cell of flag will contain a time x 1 vector of logicals (1 if
    % the frame is to be censored, 0 otherwise)
    flag = num2cell(flag,1);
    clear flag2
    clear flag3
    
    % 1 x n_subject vector containing the percentage of scrubbed frames (throughout the whole scan) 
    p_scrubbed = cell2mat(cellfun(@(x) sum(x)/length(x)*100, flag,'un',0));
    
    % If we wanted an average seed, we can either be in the
    % subject-specific seed case (we know there is only one seed type), or
    % in the case with several different seeds (up to three) to average
    % together
    if strcmp(SeedType,'Average')
        
        % According to the two options, we compute the seed time courses
        % for each subject.
        if ~is_SS
            S = CAP_get_seed_traces(tcvox,logical(sum(seed,2)),SignMat(1,:),is_SS);
        else
            S = CAP_get_seed_traces(tcvox,seed,SignMat(1,:),is_SS);
        end

        % Selects the indices appropriately (threshold or percentage
        % selection cases)
        if strcmp(SelMode,'Threshold')

            xindp = cellfun(@(x) x>T, S, 'un', 0);
            
        elseif strcmp(SelMode,'Percentage')

            T_per_act = cellfun(@(x) prctile(x,100-T), S, 'un', 0);
            xindp = cellfun(@(x,y) x>y, S,T_per_act, 'un', 0);
            
        else
            errordlg('Should never happen...');
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
        idx_sep_seeds(:,:,1) = cell2mat(xindp);
        
    % If Union or Intersection was wanted instead, then computations for
    % different seeds need to be carried independently
    else
        for idx = 1:size(seed,2)
            
            S = CAP_get_seed_traces(tcvox,logical(seed(:,idx)),SignMat(idx,:),is_SS);
            
            if strcmp(SelMode,'Threshold')

                % Computes the indexes at which we have a seed activity of interest
                xindp = cellfun(@(x) x>T, S, 'un', 0);

            elseif strcmp(SelMode,'Percentage')

                % Computes the threshold that corresponds to P percent frames for
                % each subject
                T_per_act = cellfun(@(x) prctile(x,100-T), S, 'un', 0);

                % And then uses this to select frames
                xindp = cellfun(@(x,y) x>y, S,T_per_act, 'un', 0);
                
            else
                errordlg('Should never happen...');
            end

            % flag now contains the traces with high activity AND high motion
            flag_active = cellfun(@(x,y) x & y,xindp, flag,'un',0);

            % Vector (1xn_subj) with the percentage of traces removed because of too high
            % motion and being selected as active
            p_scrubactive = cell2mat(cellfun(@(x) sum(x)/length(x)*100, flag_active,'un',0));
            
            % My indices of active/inactive frames now contain only the non
            % corrupted frames
            xindp = cellfun(@(x,y) x & ~y, xindp,flag_active,'un',0);
            
            % Updates the indices of the frames to keep
            switch SeedType
                case 'Union'
                    
                    Indices.kept.active = (Indices.kept.active) & cell2mat(xindp);
                    
                case 'Intersection'
                    
                    Indices.kept.active = (Indices.kept.active) | cell2mat(xindp);
            end
            
            idx_sep_seeds(:,:,idx) = cell2mat(xindp);
            
        end
    end
    
    % Each cell contains the frames selected as active or as inactive (if
    % inactive, the sign is reversed, i.e. inactivation is a positive
    % signal). Size: masked voxels x n_retained_frames
    Xonp = cellfun(@(x,y) x(y,:)',tcvox,num2cell(Indices.kept.active,1),'un',0);
    Xonp_full = cellfun(@(x,y) x(y,:)',tcvox,num2cell(logical(Indices.kept.active+Indices.scrubbed),1),'un',0);
    
    % Percentage of active and inactive frames retained per subject
    p_active = cell2mat(cellfun(@(x) size(x,2)/size(FDall,1)*100, Xonp,'un',0));
    
    % Matrix containing all the interesting probabilities
    p = [p_scrubbed; p_scrubactive; p_active];
end