function [CDF,Lorena] = ComputeClusteringQuality(Consensus,K_range)

    % Number of K values to check
    K = size(Consensus,3);

    % Number of pairs of frames
    n_items = size(Consensus,1);
    
    % Creates the CDF range
    c = 0:0.005:1;
    
    % Quality criterion computations are run for each explored K value...
    for k = 1:K
        
        % Sorted consensus entries
        Cons_val = sort(jUpperTriMatToVec(squeeze(Consensus(:,:,k))),'ascend');
        
        % Computation of CDF
        for i = 1:length(c)
            CDF(k,i) = sum(Cons_val <= c(i))/(n_items*(n_items-1)/2);
        end
        
        
        idx_dp = 1;
        
        for delta_perc = 0.5:0.5:10
        
            Lorena(k,idx_dp) = prctile(CDF(k,:),100-delta_perc) - prctile(CDF(k,:),delta_perc);
            idx_dp = idx_dp + 1;
        end
        
        
            
            
        % Computation of the AUC
%         AUC(k) = 0;
%         
%         for i = 2:(n_items*(n_items-1)/2)
%             AUC(k) = AUC(k) + (Cons_val(i)-Cons_val(i-1))* interp1q(c',CDF(k,:)',Cons_val(i));
%         end

        clear Cons_val
    end

%     for k = 2:K
%         
%         % Computation of Delta
%         
%         tmp_max_AUC = max(AUC(1:k-1));
% 
%         Delta(k) = (AUC(k) - tmp_max_AUC)/tmp_max_AUC;
%  
%     end
%     
%     Delta(1) = AUC(1);
    
end