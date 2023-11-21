function [hist, meanShift, posShift, negShift] = createHistogramField(Model, Field, options)
% Create the histogram of filed perturbation in each compartment   
    if ~exist('options', 'var')
        options.null = 1;
    end

    if (~isfield(options, 'create_figure') || options.create_figure == 1)
        h = figure('Name', 'Frequency histogram');
        hold on
    end
        
    if ~isfield(options, 'edges')
        options.edges = (-30 :5: 50);
    end
    
    if ~isfield(options, 'line_style')
        options.line_style = '-';
    end
    
    if ~isfield(options, 'LineWidth')
        options.LineWidth = 1;
    end
    
    if ~isfield(options, 'plot')
        options.plot = 1;
    end
    
    if ~isfield(options, 'FontSize')
        options.FontSize = 25;
    end
     if ~isfield(options, 'ylim')
        options.ylim = [0 0.3];
    end
    
    if  isfield(options, 'mask')
        Model(options.mask == 0) = -1;
    else
        options.mask = ones(size(Model));
    end
    
    listField = Field(:);
    
    nb_pixels = sum(options.mask(:));

    [hist.intra_axonal, ~] = histcounts(listField(Model == 3), options.edges);
    hist.intra_axonal = hist.intra_axonal/nb_pixels;
    
    [hist.myelin, ~] = histcounts(listField(Model == 2), options.edges);
    hist.myelin = hist.myelin/nb_pixels;
    
    [hist.extra_axonal, ~] = histcounts(listField(Model == 1), options.edges);
    hist.extra_axonal = hist.extra_axonal/nb_pixels;
    
    if  options.plot
        plot(options.edges(1:end-1), hist.extra_axonal, [options.line_style 'k'], ...
            'LineWidth', options.LineWidth)
        
        plot(options.edges(1:end-1),hist.myelin , [options.line_style 'b'], ...
            'LineWidth', options.LineWidth)
                
        plot(options.edges(1:end-1), hist.intra_axonal , [options.line_style 'r'], ...
            'LineWidth', options.LineWidth)
        
%         leg = legend('intra axonal', 'myelin', 'extra axonal');
        leg = legend('extra axonal', 'myelin', 'intra axonal');
        set(leg,'FontSize',20);
        set(leg,'fontname','Times New Roman');

%         leg = legend('intra', 'myelin', 'extra');
        xlabel('{\delta}(B) (Hz)','FontSize',25,'FontWeight','Bold','fontname','Times New Roman')
        ylabel('Frequency distribution','FontSize',25,'FontWeight','Bold','fontname','Times New Roman')
%         title('Frequency histogram')
    end
    
    if isfield(options, 'xlim')
        xlim(options.xlim);
    end
    if isfield(options, 'ylim')
        ylim(options.ylim);
    end
    
    if isfield(options, 'FontSize')
        set(gca, 'FontSize', options.FontSize);
    end
    if isfield(options, 'fontWeight')
        set(gca, 'FontWeight', options.fontWeight);
    end
     if isfield(options, 'LineWidth')
        set(gca, 'LineWidth', options.LineWidth);
    end
    
    meanShift = mean(listField(Model == 2));
    maxShift = max(listField(Model == 2));
    minShift = min(listField(Model == 2));
    
    posShift = maxShift - meanShift;
    negShift = minShift - meanShift;

end