%% Runs consensus clustering on data

% Samples the data of interest
X = SAVED.TPSelData.Act{1};

% Computes the consensus results
K_range = 2:12;
[Consensus_ordered] = CAP_ConsensusClustering(X,K_range,'subjects',0.8,100,'sqeuclidean');

[CDF,Delta] = ComputeClusteringQuality(Consensus_ordered,K_range);
