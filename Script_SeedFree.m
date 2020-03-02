%% This is an example script to run CAP analysis without the GUI
% In this script, we assume that the user is interested in performing a
% seed-free analysis with one population of subjects



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



%% 2. Specifying the main parameters

% Threshold of FD above which to scrub out the frame and also the t-1 and
% t+1 frames (if you want another scrubbing setting, directly edit the
% code)
Tmot = 0.5;

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
% Note that since this is a seed-free analysis, we do not need any
% seed-related information to be given as argument
[Xon,~,Indices] = CAP_find_activity_SeedFree(TC,FD,Tmot); 



%% 4. Consensus clustering (if wished to determine the optimum K)

% This specifies the range of values over which to perform consensus
% clustering: if you want to run parallel consensus clustering processes,
% you should feed in different ranges to each call of the function
K_range = 2:50;

% Have each of these run in a separate process on the server =)
[Consensus] = CAP_ConsensusClustering(Xon,K_range,'items',Pcc/100,N,'correlation');

% Calculates the quality metrics
[~,Qual] = ComputeClusteringQuality(Consensus,[]);

% Qual should be inspected to determine the best cluster number(s)

% You should fill this with the actual value 
K_opt = 



%% 5. Clustering into CAPs

[CAP,~,~,idx] = Run_Clustering(cell2mat(Xon),...
        K_opt,mask,brain_info,Pp,Pn,n_rep,[],SeedType);

    

%% 6. Computing metrics

% The TR of your data in seconds
TR = 

[ExpressionMap,Counts,Entries,Avg_Duration,Duration,TransitionProbabilities,...
    From_Baseline,To_Baseline,Baseline_resilience,Resilience,Betweenness,...
    InDegree,OutDegree,SubjectEntries] = Compute_Metrics_simpler(idx,...
    Indices.kept.active,Indices.scrubbedandactive,K_opt,TR);