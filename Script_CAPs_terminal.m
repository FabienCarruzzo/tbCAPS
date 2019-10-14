%% This is an example script to run the CAPs analyses without the GUI
% Coded by Thomas for Elvira



%% 1. Loading the data files

% Data: cell array, each cell of size n_TP x n_masked_voxels
TC =

% Mask: n_voxels x 1 logical vector
mask =

% Header: the header (obtained by spm_vol) of one NIFTI file with proper
% data dimension and .mat information
brain_info =

% Framewise displacement: a n_TP x n_subj matrix with framewise
% displacement information
FD =

% Seed: a n_masked_voxels x n_seed logical vector with seed information
Seed = 


%% 2. Specifying the main parameters

% Threshold above which to select frames
T = 1;

% Selection mode ('Threshold' or 'Percentage')
SelMode = 'Threshold';

% Threshold of FD above which to scrub out the frame and also the t-1 and
% t+1 frames (if you want another scrubbing setting, directly edit the
% code)
Tmot = 0.3;

% Type of used seed information: select between 'Average','Union' or
% 'Intersection'
SeedType = 'Intersection';

% Contains the information, for each seed (each row), about whether to
% retain activation (1 0) or deactivation (0 1) time points
SignMatrix = [0 1; 1 0; 1 0];

isSeedSubjectSpecific = 0;

% Percentage of positive-valued voxels to retain for clustering
Pp = 100;

% Percentage of negative-valued voxels to retain for clustering
Pn = 100;

% Number of repetitions of the algorithm
n_rep = 50;


%% 3. Selecting the frames to analyse

[Xonp,p,Indices,idx_sep_seed] = CAP_find_activity(Data_train,...
    SAL,T,FD_train,Tmot,SelMode,SeedType,SignMatrix,isSeedSubjectSpecific);
    
    
%% 4. Consensus clustering (if wished to determine the optimum K)

% Computes the consensus results
K_range = 2:50;
[Consensus] = CAP_ConsensusClustering(Xonp,K_range,'items',0.8,20,'sqeuclidean');

% Calculates the quality metrics
[~,QC] = ComputeClusteringQuality(Consensus,K_range);

% Here, you should plot QC and check what value of K yields the
% best value
figure; bar(K_range,QC);

K_opt = 


%% 5. Clustering into CAPs

[CAP,~,STDCAP,idx,CorrDist,sfrac] = Run_Clustering(cell2mat(Xonp),...
    K_opt,mask,brain_info,Pp,Pn,n_rep,idx_sep_seed,SeedType);


%% 6. Computing metrics

TR = 

[TPM,Counts,Number,Avg_Duration,Duration,TM,From_Baseline,...
    To_Baseline,Baseline_resilience,Resilience,Betweenness,kin,kout] = ...
    Compute_Metrics_simpler(idx,Indices.kept.active,...
    Indices.scrubbedandactive,K_opt,TR);