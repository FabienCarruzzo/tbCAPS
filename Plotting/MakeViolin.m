%% Makes a violin plot with overlapped boxplot, tuning colors
% Y and O are vectors containing the data to plot for the two conditions
% Color has light colors for background (first cell) and dark colors for
% average (second cell), for maximum 4 populations
function [box,h_viol,ah] = MakeViolin(Y,ah,Lab,YLabel,Color,n_pop,n_states)
    
    % Range of the plot
    Max = max(max((Y)));
    Min = min(min((Y)));
    Max = Max + 0.15*max(abs(Min),Max);
    Min = Min - 0.15*max(abs(Min),Max);
              
    % Plots the distribution
    ylim(ah,[Min,Max]);

    % Colors for the background
    Col_final = num2cell(Color{1},2);
    Col_final = Col_final(1:n_pop);
    Col_final = repmat(Col_final,n_states,1);
    
    % Colors for the foreground
    Col_final2 = num2cell(Color{2},2);
    Col_final2 = Col_final2(1:n_pop);
    Col_final2 = repmat(Col_final2,n_states,1);
    
    Lab_final = {};
    
    for i = 1:n_states
        Lab_final = [Lab_final, repmat({Lab{i}},1,n_pop)];
    end
    
    h_viol = distributionPlot(ah,Y','showMM',0,'color',Col_final,'yLabel',YLabel,'xNames',Lab_final);
    h_viol=0;
    hold(ah,'on');
    box = boxplot(Y','Parent',ah,'Labels',Lab_final);
    
    h_box = findobj(box,'Tag','Box');
    set(h_box,'color','k','LineWidth',2);
    
    h_median = findobj(box,'Tag','Median');
    h_outliers = findobj(box,'Tag','Outliers');
    
    for i = 1:n_pop
        for j = 1:n_states
    
            idx_oi = i+(j-1)*n_pop;
            
            set(h_median(idx_oi),'color',Col_final2{idx_oi},'LineWidth',2);
            set(h_outliers(idx_oi),'Marker','o','MarkerFaceColor',Col_final2{n_pop*n_states-idx_oi+1},...
                'MarkerEdgeColor',Col_final2{idx_oi},'MarkerSize',3);

        end
    end
    
    set(ah,'Box','off');
end