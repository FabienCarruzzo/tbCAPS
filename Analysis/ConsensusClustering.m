%% This function performs consensus clustering over a range of K values
% The goal is to provide a measure of how good each value of K is
% 
% Inputs:
% - X is the data matrix (n_DP x n_DIM)
% - K_range is the range of K values to examine
% - Subsample_type defines how subsampling is done: across items (data
% points) if 'items', and across dimensions if 'dimensions'
% - Subsample_fraction is the fraction of the original data points, or
% dimensions, to keep for a given fold
% - n_folds is the number of folds over which to run
function [Consensus_ordered] = ConsensusClustering(X,K_range,Subsample_type,Subsample_fraction,n_folds,DistType)

    % Number of data points
    n_items = size(X,1);
    
    % Number of dimensions
    n_dims = size(X,2);
    
    Consensus = zeros(n_items,n_items,length(K_range));
    Consensus_ordered = zeros(n_items,n_items,length(K_range));
    
    % Loop over all K values to assess
    for k = 1:length(K_range)
    
        disp(['Running consensus clustering for K = ',num2str(K_range(k)),'...']);
        
        % Connectivity matrix that will contain 0s or 1s depending on whether
        % elements are clustered together or not
        M = zeros(n_items,n_items,n_folds);
        I = zeros(n_items,n_items,n_folds);
        
        disp('before h loop');
        
        % Loops over all the folds to perform clustering for
        for h = 1:n_folds
            
            switch Subsample_type
                case 'items'
                    
                    % Number of items to subsample
                    n_items_ss = floor(Subsample_fraction*n_items);
                    
                    % Does the subsampling
                    [X_ss,tmp_ss] = datasample(X,n_items_ss,1,'Replace',false);
                    
                    % Vector
                    I_vec = zeros(n_items,1);
                    I_vec(tmp_ss) = 1;
                    
                    % Constructs the indicator matrix
                    for i = 1:length(I_vec)
                        for j = 1:length(I_vec)
                            if (I_vec(i) == I_vec(j)) && (I_vec(i) > 0)
                                I(i,j,h) = 1;
                            end
                        end
                    end
                    
                case 'dims'
                    
                    % Number of dimensions to subsample
                    n_dims_ss = floor(Subsample_fraction*n_dims);
                    
                    % Does the subsampling
                    [X_ss,tmp_ss] = datasample(X,n_dims_ss,2,'Replace',false);
                    
                    % Constructs the indicator matrix
                    I(:,:,h) = ones(n_items,n_items);
                    
                otherwise
                    errordlg('PROBLEM IN TYPE OF SUBSAMPLING');
            end
            
            % Does the clustering (for now, only with k-means), so that IDX
            % contains the indices for each datapoint
            IDX = kmeans(X_ss,K_range(k),'Distance',DistType,'Replicates',1,'Start','uniform');
            
            % Builds the connectivity matrix
            M(:,:,h) = Build_Connectivity_Matrix(IDX,tmp_ss,Subsample_type,n_items); 
            
            clear I_vec
            clear X_ss
            clear tmp_ss
            clear IDX
        end
        
        % Constructs the consensus matrix for the considered K
        Consensus(:,:,k) = sum(M,3)./sum(I,3); 
        
        tree = linkage(squeeze(1-Consensus(:,:,k)),'average');

        % Leaf ordering to create a nicely looking matrix
        leafOrder = optimalleaforder(tree,squeeze(1-Consensus(:,:,k)));
        
        % Ordered consensus matrix
        Consensus_ordered(:,:,k) = Consensus(leafOrder,leafOrder,k);
        
        clear leafOrder
        clear Dist_vec
        clear test
        clear IDX
        clear M
        clear I
    end
end