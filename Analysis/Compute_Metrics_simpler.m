%% This function computes all the metrics of the CAP analysis framework
% Inputs:
%
% - idx depicts the indices of the frames that have been clustered (the
% cluster to which they belong)
% - xindp1 depicts the indices of the frames that have been selected as
% activation moments
% - xindn1 is the same for deactivation time points
% - sind depicts scrubbed frames
% - n_clusters is the number of clusters used for state disentanglement
% - TR is the TR of the experiment
% - CAPType denotes the type of clustering that has been chosen (only
% activation, only deactivation, or both)
% 
% Outputs:
% 
% - TPM (n_subjects x n_frames) is the state sequence matrix
% - Counts (-> raw/frac -> scrubbed/baseline/state: n_subj x n_states)
% contains the counts (raw and normalized)
% - Number (n_subj x n_seqtype) with sequence type: scrub, act/deact/-,
% baseline, then states
% - Avg_Duration (n_subj x n_seqtype) contains the average duration of a
% state for a subject
% - Duration (n_subj array of size 1 x n_sequences) contains the duration
% of the state sequences for each subject
% - TM (n_seqtype x n_seqtype) encompasses the transitions that exist
function [TPM,Counts,Number,Avg_Duration,Duration,TM,From_Baseline,...
    To_Baseline,Baseline_resilience,Resilience,Betweenness,kin,kout,SubjectEntries] = ...
    Compute_Metrics_simpler(idx,xindp1,sind,n_clusters,TR)

    %% Declaration of parameters

    % Number of subjects
    n = size(xindp1,2);
    
    % Number of frames (full)
    n_frames = size(xindp1,1);
    
    % Cumulative frame counts across subjects (1xn_subjects) for frames
    dd2p = cumsum([1,sum(xindp1,1)]);
    
    % TPM will contain the sequence of states for all subjects
    TPM = zeros(n,n_frames);

    % I want a vector that specifies, for each subject, how many frames
    % from each CAP are present
    SubjectEntries = zeros(n,n_clusters+3);
    
    % Filling of TP for each subject
    for i = 1:n

        % We set scrubbed time points at -1 if we want to discard scrubbed
        % frames; else, we set them to the retrieved indices (and -1
        % denotes unassigned cases)
        TPM(i,xindp1(:,i)) = idx(dd2p(i):dd2p(i+1)-1)'; 
        TPM(i,sind(:,i)) = -1;
    end
    
    for kk = -1:n_clusters+1
        SubjectEntries(:,kk+2) = sum(TPM == kk,2);
    end
    
    % Computation of the counts
    [Counts,Duration,Number,Avg_Duration,TM] = CAP_ComputeMetrics(TPM,n_clusters+1,TR,n_frames);
    
    %% Computation of novel metrics
    
    % I want to consider: all my K CAPs, and the baseline state (but not
    % the (K+1)th Unassigned CAP or scrubbed frames)
    TP = TM(3:end-1,3:end-1,:);
    
    % How often do I reach a CAP from baseline
    From_Baseline = squeeze(TM(2,3:end-1,:))';
    
    % How often do I reach the baseline from a CAP
    To_Baseline = squeeze(TM(3:end-1,2,:))';
    
    % How long do I stay in baseline
    Baseline_resilience = squeeze(TM(2,2,:));

    % First, I can compute the probability to stay within a given state; I
    % then set the probability to 0 for my next computations
    for s = 1:size(TP,3)
        for k = 1:n_clusters

            Resilience(s,k) = TP(k,k,s);
            
            TP(k,k,s) = 0;
        end

        Betweenness(s,:) = betweenness_wei(squeeze(TP(:,:,s)));
        [kin(s,:),kout(s,:)] = degrees_dir(squeeze(TP(:,:,s)));

    end
    
end