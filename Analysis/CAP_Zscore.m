function [iCAPs_z] = CAP_Zscore(iCAPs)
    for i = 1:size(iCAPs,1)
        [a,b] = hist(iCAPs(i,:), 100);
        aind = find(a == max(a));
        med  = b(aind(1));
        iCAPs_z(i,:) = (iCAPs(i,:)-med)/sqrt((sum((iCAPs(i,:)-med).^2))/length(iCAPs(i,:)));    % normalization copied from Isik
    end
end