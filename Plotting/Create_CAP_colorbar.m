%% Creates a nice looking colorbar for the CAP display
function [Ah,h] = Create_CAP_colorbar(absmin,absmax,absstep,colT,lab,Ah,Ori,CB1,CB2,n_steps)

    H_range = [absmin absmax]; % The colormap is symmetric around zero

    % Set the Min/Max T-values for alpha coding
    A_range = [0 1];
    
    % Set the labels for the colorbar
    hue_label = lab;
    
    colrange = linspace(absmin,absmax,256);
    
    switch Ori
        case 'Horizontal'
            
            y = linspace(A_range(1), A_range(2), 256); 
            % x represents the range in alpha (abs(t-stats))
            x = linspace(H_range(1), H_range(2), 256);
            % y represents the range in hue (beta weight difference)
            [X,Y] = meshgrid(x,y); % Transform into a 2D matrix

            h=imagesc(x,y,X,'Parent',Ah); 
            axis(Ah,'xy'); % Plot the colorbar
            set(Ah, 'Xcolor', 'k', 'Ycolor', 'k','YTickLabel','','YTick',[],'XTick',absmin:absstep:absmax,'FontSize',8);
            set(Ah, 'XAxisLocation', 'bottom');
            xlabel(Ah,hue_label,'FontSize',8);

            A = ones(size(X));
            A(abs(X) < colT) = 0;
            A = reshape(A,256,256);
            
        case 'Vertical'
            
            x = linspace(A_range(1), A_range(2), 256); 
            % x represents the range in alpha (abs(t-stats))
            y = linspace(H_range(1), H_range(2), 256);
            % y represents the range in hue (beta weight difference)
            [X,Y] = meshgrid(x,y); % Transform into a 2D matrix

            h=imagesc(x,y,Y,'Parent',Ah); 
            axis(Ah,'xy'); % Plot the colorbar
            set(Ah, 'Xcolor', 'k', 'Ycolor', 'k','XTickLabel','','XTick',[],'YTick',absmin:absstep:absmax,'FontSize',8);
            set(Ah, 'YAxisLocation', 'right');
            ylabel(Ah,hue_label,'FontSize',8);

            A = ones(size(Y));
            A(abs(Y) < colT) = 0;
            A = reshape(A,256,256);
    end
    
        tmp_cmap = cbrewer(CB1,CB2,n_steps);
        colormap(Ah,flipud(tmp_cmap));
    
     set(h,'AlphaData',A);
end