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
function [Metrics,Metrics_ext] = ...
    Compute_Metrics(idx,xindp1,sind,n_clusters,TR,CAP,X_extended,T_assign,d)

    %% Declaration of parameters

    % Number of subjects
    n = size(xindp1,2);
    
    % Number of frames (full)
    n_frames = size(xindp1,1);
    
    % Cumulative frame counts across subjects (1xn_subjects) for frames
    dd2p = cumsum([1,sum(xindp1,1)]);
    
    % TPM will contain the sequence of states for all subjects
    TPM = zeros(n,n_frames);
    
    % Sequence of states if we try to also match the scrubbed frames to the
    % CAPs
    TPM_extended = zeros(n,n_frames);

    
    % Filling of TP for each subject
    for i = 1:n
        
        % Indices of the CAPs to which scrubbed frames would belong for the
        % subject of interest (in is, a vector)
        is = CAP_AssignFrames(CAP,X_extended{i},d,T_assign);

        % We set scrubbed time points at -1 if we want to discard scrubbed
        % frames; else, we set them to the retrieved indices (and -1
        % denotes unassigned cases)
        TPM(i,xindp1(:,i)) = idx(dd2p(i):dd2p(i+1)-1)'; 
        TPM_extended(i,xindp1(:,i)) = idx(dd2p(i):dd2p(i+1)-1)'; 

        TPM(i,sind(:,i)) = -1;
        
        if ~isempty(is)
            is(is > n_clusters) = -1;
            TPM_extended(i,sind(:,i)) = is;
        end
   
    end
    
    % Computation of the counts
    [Counts,Duration,Number,Avg_Duration,TM] = CAP_ComputeMetrics(TPM,n_clusters+1,TR,n_frames);
    [Counts_ext,Duration_ext,Number_ext,Avg_Duration_ext,TM_ext] = CAP_ComputeMetrics(TPM_extended,n_clusters+1,TR,n_frames);
    
    %% Computation of novel metrics
    
    % I want to consider: all my K CAPs, and the baseline state
    TP = TM(3:end-1,3:end-1,:);
    TP_ext = TM_ext(3:end-1,3:end-1,:);
    
    % How often do I reach a CAP from baseline
    From_Baseline = squeeze(TM(2,3:end-1,:))';
    From_Baseline_ext = squeeze(TM_ext(2,3:end-1,:))';
    
    % How often do Ireach the baseline from a CAP
    To_Baseline = squeeze(TM(3:end-1,2,:))';
    To_Baseline_ext = squeeze(TM_ext(3:end-1,2,:))';
    
    % How long do I stay in baseline
    Baseline_resilience = squeeze(TM(2,2,:));
    Baseline_resilience_ext = squeeze(TM_ext(2,2,:));

    % First, I can compute the probability to stay within a given state; I
    % then set the probability to 0 for my next computations
    for s = 1:size(TP,3)
        for k = 1:n_clusters

            Resilience(s,k) = TP(k,k,s);
            Resilience_ext(s,k) = TP_ext(k,k,s);
            
            TP(k,k,s) = 0;
            TP_ext(k,k,s) = 0;
        end

        Betweenness(s,:) = betweenness_wei(squeeze(TP(:,:,s)));
        [kin(s,:),kout(s,:)] = degrees_dir(squeeze(TP(:,:,s)));

        Betweenness_ext(s,:) = betweenness_wei(squeeze(TP_ext(:,:,s)));
        [kin_ext(s,:),kout_ext(s,:)] = degrees_dir(squeeze(TP_ext(:,:,s)));
    end

    Metrics.From_Baseline = From_Baseline;
    Metrics.To_Baseline = To_Baseline;
    Metrics.Baseline_resilience = Baseline_resilience;
    Metrics.Resilience = Resilience;
    Metrics.TP = TP;
    Metrics.Betweenness = Betweenness;
    Metrics.kin = kin;
    Metrics.kout = kout;
    Metrics.TPM = TPM;
    Metrics.Counts = Counts;
    Metrics.Duration = Duration;
    Metrics.Avg_Duration = Avg_Duration;
    Metrics.Number = Number;
    
    Metrics_ext.From_Baseline = From_Baseline_ext;
    Metrics_ext.To_Baseline = To_Baseline_ext;
    Metrics_ext.Baseline_resilience = Baseline_resilience_ext;
    Metrics_ext.Resilience = Resilience_ext;
    Metrics_ext.TP = TP_ext;
    Metrics_ext.Betweenness = Betweenness_ext;
    Metrics_ext.kin = kin_ext;
    Metrics_ext.kout = kout_ext;
    Metrics_ext.TPM = TPM_extended;
    Metrics_ext.Counts = Counts_ext;
    Metrics_ext.Duration = Duration_ext;
    Metrics_ext.Avg_Duration = Avg_Duration_ext;
    Metrics_ext.Number = Number_ext;
    
end