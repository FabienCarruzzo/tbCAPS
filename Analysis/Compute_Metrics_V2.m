function [TPM,TM] = Compute_Metrics_V2(idx,xindp1,sind,n_clusters,TR)

    % Number of subjects
    n = size(xindp1,2);
    
    % Number of frames
    n_frames = size(xindp1,1);
    
    % This will contain our transitions
    TM = zeros(n_clusters+2,n_clusters+2,n);
    
    % Cumulative frame counts across subjects (1xn_subjects) for frames
    dd2p = cumsum([1,sum(xindp1,1)]);
    
    % TPM will contain the sequence of states for all subjects
    TPM = zeros(n,n_frames);
    
    % Computation of the state matrix
    
    % Filling of TP for each subject
    for i = 1:n
        % We set scrubbed time points at -1
        TPM(i,xindp1(:,i)) = idx(dd2p(i):dd2p(i+1)-1)'; 
        TPM(i,sind(:,i)) = -1;
    end
    
    % We now want to compute dynamical metrics for each subject
    for j = 1:n

        % For each state transition, we increment properly the matrix of
        % transitions
        for k = 1:n_frames-1  
            TM(TPM(j,k)+2,TPM(j,k+1)+2,j) = TM(TPM(j,k)+2,TPM(j,k+1)+2,j) + 1;
        end
    end
    
    % We normalize the matrix by all the transitions
    TM = TM/(n_frames-1);
    
    
end