%% This function constructs a connectivity matrix from the indices of
% cluster assignment
function [M] = Build_Connectivity_Matrix(IDX,tmp_ss,type,n_items)

    IDX_full = zeros(n_items,1);

    switch type
        case 'items'
            IDX_full(tmp_ss) = IDX;
        case 'dims'
            IDX_full = IDX;
        case 'subjects'
            IDX_full(tmp_ss) = IDX;
    end
    
    M = zeros(n_items,n_items);

    for i = 1:length(IDX_full)
        for j = 1:length(IDX_full)
            if (IDX_full(i) == IDX_full(j)) && (IDX_full(i) > 0)
                M(i,j) = 1;
            end
        end
    end
end