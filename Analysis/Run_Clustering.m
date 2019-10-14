function [CP2,Disp,Std_Clusters,idx,d,sfrac] = Run_Clustering(XONn,n_clusters,mask,brain_info,maskP,maskN,n_rep,idx_sep_seeds,SeedType)

    % Number of seeds
    n_seeds = size(idx_sep_seeds,3);
    
    % Number of subjects
    n_subjects = size(idx_sep_seeds,2);

    % Number of possible seed combinations
    switch n_seeds
        case 1
            n_combos = 1;
            combvec = [1];
        case 2
            n_combos = 3;
            combvec = [1 0; 0 1; 1 1];
        case 3
            n_combos = 7;
            combvec = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1];
        otherwise
            errordlg('PROBLEM AT SEED FRACTIONS');
    end

    % Will contain the fraction of frames linked to a given seed (if using
    % the 'Intersection' method)
    sfrac = zeros(n_subjects,n_clusters,n_combos);

    % 'Filtering so that we only use the largest activation and deactivation
    % spots for clustering
    XONn_filtered = CAP_mask4kmeans(XONn,maskP,maskN,6,mask,brain_info);
       
    % Rows datapoints, columns variable (so here every 70700 activation is
    % a datapoint)
    % idx will contain 1462 elements (the index of the cluster to which the
    % considered datapoint belongs
    [idx,CP] = kmeans(XONn_filtered',n_clusters,'distance','correlation','replicates',n_rep,'empty','drop','maxiter',100,'Display','iter');
        
    % idx2counts is of size K (number of clusters) and has the number of
    % datapoints classified within a given cluster)

    % disp('idx2counts:');
    idx2counts = histc(idx, 1:max(idx));

    % Output = Input(IX)
    [~,IX] = sort(idx2counts,'descend');

    % Size Kx70700 (location of each cluster); clusters are put with 'the
    % most prominent one first'
    CP = CP(IX,:); idx2 = idx; % order by occurrence

    % Changes the datapoint indexes so that they fit the new clusters order
    for l=1:max(idx), idx2(idx==IX(l))=l; end
    idx=idx2;

    CP2 = zeros(n_clusters,size(CP,2));
    Disp = zeros(1,n_clusters);
    Std_Clusters = zeros(size(CP,2),n_clusters);

    % For each cluster index
    for l=1:max(idx)       
        % Averages all data points belonging to one specific cluster, and
        % stores the obtained pattern as a cell in CP2
        CP2(l,:) = mean(XONn(:,idx==l),2); %./ ( std(XON(:,idx==l),[],2) ) * sqrt(length(idx==l)); % Liu&Duyn

        % Measure of dispersion within the cluster considered
        Disp(l) = mean(corr(CP2(l,:)',XONn(:,idx==l)));

        Std_Clusters(:,l) = std(XONn(:,idx==l),[],2);
    end
    
    % d contains the correlation values of all frames to the CAPs
    r = corr(CP2',XONn);

    d = zeros(1,length(idx));
    for k=1:max(idx)
        d(idx==k) = r(k,idx==k);
    end
    
    % Added part to compute the fraction of frames assigned to a given seed
    % if using the intersection option (in which a data point is retained
    % as long as at least one seed region becomes significantly (de)active)
    if strcmp(SeedType,'Intersection')
        
        % idx_all will contain the clustering indices as put on a whole
        % temporal scale (time x subjects)
        idx_all = zeros(size(idx_sep_seeds,1),n_subjects);
        
        % Index to properly fill in the matrix by adding up the number of
        % frames per subject every time
        tmp_loc = 1;
        

        for s = 1:n_subjects
            tmp = sum(squeeze(idx_sep_seeds(:,s,:)),2);
            tmp(tmp >= 1) = 1;
            tmp = logical(tmp);
            idx_all(tmp,s) = idx(tmp_loc:(tmp_loc+sum(tmp)-1));
            tmp_loc = tmp_loc + sum(tmp);
        end
            
        
        
        % I will compute my fractions for each seed combination of
        % interest; for example, 3 seeds would yield 7 possible
        % combinations
        for s = 1:n_subjects
            for t = 1:size(idx_sep_seeds,1)
                
                tmp = squeeze(idx_sep_seeds(t,s,:));   
                
                % If there is at least one seed active or deactive at
                % the point of interest, we update the sfrac count at the
                % appropriate CAP
                if sum(tmp) > 0
                    sfrac(s,idx_all(t,s),find(ismember(combvec,tmp','rows'))) = ...
                        sfrac(s,idx_all(t,s),find(ismember(combvec,tmp','rows'))) + 1;   
                end
            end
        end 
    end 
end