%% This function checks if the data that have been entered are of consistant
% dimensions
function [IsOK,Problems] = CAP_IsDataOK(TC,FD,mask,BI)

    % By default, stuff is not OK
    IsOK = 0;

    % If the dimensions of TC are consistant...
    if size(TC,1) == 1 && size(TC,2) >...
                3 && sum(cell2mat(cellfun(@(x) size(x,1) == size(TC{1},1) && ...
                size(x,2) == size(TC{1},2),TC,'un',0))) == size(TC,2)

        % If the dimensions of FD fit the dimensions of TC...
        if size(FD,1) == size(TC{1},1) && size(FD,2) == size(TC,2)

            % If the dimensions of the mask fit the rest...
            if size(mask,2) == 1 && sum(mask,1) == size(TC{1},2)

                % If the dimensions of the brain info are consistant...
                if isfield(BI,'dim') && isfield(BI,'mat') && size(BI.mat,1) == 4 && size(BI.mat,2) == 4 && size(BI.dim,1) == 1 && size(BI.dim,2) == 3

                    % Then, everything is OK !
                    IsOK = 1;
                    Problems = 'No problem !!';
                else
                    Problems = 'Inconsistant brain information dimensions';
                end
            else
                Problems = 'Inconsistant mask dimensions compared to time courses';
            end
        else
            Problems = 'Inconsistant dimensions between time courses and motion file';
        end
        
    else
        Problems = 'Inconsistant time courses dimensions';
    end


end