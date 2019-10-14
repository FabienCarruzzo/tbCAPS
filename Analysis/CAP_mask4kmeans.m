%% This function sets to zero the value of the voxels that are judged noise
% X is the data to mask (n_vox x n_frames)
% topP and bottomP are the percentages of data retained (top and lowest
% activity ones respectively) by the masking procedure
% n_cv is the number of voxels in a group that must be reached for that
% group not to be masked
function [Y] = CAP_mask4kmeans(X,topP,bottomP,n_cv,mask,ai)
    
    % We want to exit with a matrix of the same size as X, for which each
    % frame has had its noise-related elements set to zero
    Y = zeros(size(X));

    % For all frames...
    for i = 1:size(X,2)
        % We sort in descending and ascending order to get the indexes of
        % the considered frame matching the top percentage of high or low
        % activity
        [~,Isorted] = sort(X(:,i),'descend');
        [~,Isortedrev] = sort(X(:,i),'ascend');
        
        % Contains the indexes of all the points with activity of interest,
        % both high and low
        I = [Isorted(1:round(topP/100*length(Isorted)))',Isortedrev(1:round(bottomP/100*length(Isortedrev)))'];
        
        % If the elements of the frame belong to the indexes, X_binary has
        % the corresponding entry set to 1. Else, it is set to 0
        X_binary = ismember(1:length(X(:,i)),I);
        
        % To perform the opening operation, we must convert X_binary to a
        % 3D volume
        temp = nan(size(mask)); 
        temp(mask) = X_binary; 
        temp(isnan(temp)) = 0;
        temp = reshape(temp,ai.dim);
        
        % At this state, temp is a binary image on which the small clusters
        % of neigboring voxels retained (less than n_cv neighbours) have
        % been removed
        temp = bwareaopen(temp,n_cv);
        
        % We convert temp back to a 1D vector, and conserve only the
        % elements of interest (70700 voxels) with 0 and 1 respectively
        % denoting 'not to consider for k-means' and 'to consider for
        % k-means'
        temp = temp(:);
        Y(:,i) = temp(mask);
    end
    
    % Y is a binary matrix of same size as X; multiplying element by
    % element, we set to zero the elements of X that we want to neglect for
    % clustering
    Y = Y.*X;
end