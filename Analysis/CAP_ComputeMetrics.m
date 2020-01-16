function [Counts,Duration,Number,Avg_Duration,TM] = CAP_ComputeMetrics(TPM,n_clusters,TR,n_frames)

    TM = zeros(n_clusters+2,n_clusters+2,size(TPM,1));

    % We first compute the amount and percentage of scrubbed frames
    % Raw counts
    Counts.raw.scrubbed = sum(TPM == -1,2);
    Counts.raw.notpicked = sum(TPM == 0,2);

    % Normalized counts
    Counts.frac.scrubbed = 100*sum(TPM == -1,2)./sum(TPM,2);
    Counts.frac.notpicked = 100*sum(TPM == 0,2)./sum(TPM,2);

    for idx_state = 1:n_clusters
        Counts.raw.state(:,idx_state) = sum(TPM == idx_state,2);
        Counts.frac.state(:,idx_state) = 100*sum(TPM == idx_state,2)./sum(TPM > 0,2);
    end

    % Computation of state number and duration

    % We now want to compute dynamical metrics for each subject
    for j = 1:size(TPM,1)

        % Marks the indices when we switch state (last element is idx_max+1)
        TS = find(diff([-1337 TPM(j,:) -1337])~=0);

        % Length of each period (in s.)
        Length = diff(TS)*TR;

        % Starting index of each period (in s., first element at 0 s.)
        StartTime = (TS(1:end-1)-1)*TR;

        % Type of each period
        SeqType = TPM(j,TS(1:end-1)); 

        % Number of different states that have been entered
        n_Periods = length(Length);

        % We now want to actually count how many times we enter each state,
        % and how long we stay in each state
        n_Sequences = zeros(1,n_clusters+2);
        length_Sequences = cell(1,n_clusters+2);

        % We go through all sequences
        for i = 1:n_Periods

            % We increase the appropriate counter
            n_Sequences(SeqType(i)+2) = n_Sequences(SeqType(i)+2)+1;
            length_Sequences{SeqType(i)+2} = [length_Sequences{SeqType(i)+2},Length(i)]; 
        end

        for ns = 1:n_clusters+2
            % Average duration of states
            avg_length_Sequences(ns) = mean(length_Sequences{ns});
        end

        Duration{j} = length_Sequences;
        Number(j,:) = n_Sequences;
        Avg_Duration(j,:) = avg_length_Sequences; 

        % For each state transition, we increment properly the matrix of
        % transitions
        for k = 1:n_frames-1  
            TM(TPM(j,k)+2,TPM(j,k+1)+2,j) = TM(TPM(j,k)+2,TPM(j,k+1)+2,j) + 1;
        end
    end

    % We normalize the matrix by all the transitions
    TM = TM/(n_frames-1);
    
end