function [i] = CAP_AssignFrames(Cp,XON,d,T)

    % Value of correlation below which 5% of all values lie (so gives an
    % estimate of a 'bad correlation for a trace belonging to a cluster'
    CT = prctile(d,T);

    % cluster assignment in young
    r = corr(Cp',XON);
    [c,i] = max(r);

    % If the correlation value is too low (below threshold), the index is set
    % to a new non-existing group
    i(c<CT) = size(Cp,1)+1;
end