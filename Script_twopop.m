%% This is an example script to run the CAPs analyses without the GUI
% In this script, we assume one population of subjects only



%% 1. Loading the data files

% Data: cell array, each cell of size n_TP x n_masked_voxels
TC1 =
TC2 =

% Mask: n_voxels x 1 logical vector
mask =

% Header: the header (obtained by spm_vol) of one NIFTI file with proper
% data dimension and .mat information
brain_info =

% Framewise displacement: a n_TP x n_subj matrix with framewise
% displacement information
FD1 =
FD2 =

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
Tmot = 0.5;

% Type of used seed information: select between 'Average','Union' or
% 'Intersection'
SeedType = 'Average';

% Contains the information, for each seed (each row), about whether to
% retain activation (1 0) or deactivation (0 1) time points
SignMatrix = [1,0];

% Percentage of positive-valued voxels to retain for clustering
Pp = 100;

% Percentage of negative-valued voxels to retain for clustering
Pn = 100;

% Number of repetitions of the K-means clustering algorithm
n_rep = 50;

% Percentage of frames to use in each fold of consensus clustering
Pcc = 80;

% Number of folds we run consensus clustering for
N = 50;



%% 3. Selecting the frames to analyse    

% Xon will contain the retained frames, and Indices will tag the time
% points associated to these frames, for each subject (it contains a
% subfield for retained frames and a subfield for scrubbed frames)
[Xon1,~,Indices1] = CAP_find_activity(TC1,seed,T,FD1,Tmot,SelMode,SeedType,SignMatrix);
[Xon2,~,Indices2] = CAP_find_activity(TC2,seed,T,FD2,Tmot,SelMode,SeedType,SignMatrix);
    

%% 4. Consensus clustering (if wished to determine the optimum K)

% This specifies the range of values over which to perform consensus
% clustering: if you want to run parallel consensus clustering processes,
% you should feed in different ranges to each call of the function
K_range = 2:20;

% Have each of these run in a separate process on the server =)
[Consensus] = CAP_ConsensusClustering(Xon1,K_range,'items',Pcc/100,N,'correlation');

% Calculates the quality metrics
[~,Qual] = ComputeClusteringQuality(Consensus,[]);

% Qual should be inspected to determine the best cluster number(s)

% You should fill this with the actual value 
K_opt = 



%% 5. Clustering into CAPs

[CAP,~,~,idx1,CorrDist] = Run_Clustering(cell2mat(Xon1),...
        K_opt,mask,brain_info,Pp,Pn,n_rep,[],SeedType);


    
%% 6. Assignment of the frames from population 2

% Parameter that governs the stringency of assignment: if Ap = 5%, we
% assign a frame to a CAP if spatial correlation exceeds the 5th percentile
% of the distribution of spatial correlations between the CAP, and its
% constituting frames
Ap = 5;

idx2 = CAP_AssignFrames(CAP,cell2mat(Xon2),CorrDist,Ap)';



%% 7. Computing metrics

% The TR of your data in seconds
TR = 

[ExpressionMap1,Counts1,Entries1,Avg_Duration1,Duration1,TransitionProbabilities1,...
    From_Baseline1,To_Baseline1,Baseline_resilience1,Resilience1,Betweenness1,...
    InDegree1,OutDegree1,SubjectEntries1] = Compute_Metrics_simpler(idx1,...
    Indices1.kept.active,Indices1.scrubbedandactive,K_opt,TR);

[ExpressionMap2,Counts2,Entries2,Avg_Duration2,Duration2,TransitionProbabilities2,...
    From_Baseline2,To_Baseline2,Baseline_resilience2,Resilience2,Betweenness2,...
    InDegree2,OutDegree2,SubjectEntries2] = Compute_Metrics_simpler(idx2,...
    Indices2.kept.active,Indices2.scrubbedandactive,K_opt,TR);