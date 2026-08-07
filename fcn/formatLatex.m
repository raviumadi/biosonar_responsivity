function style = formatLatex(fig, preset)
%FORMATLATEX Apply manuscript-ready sizing and typography to a figure.
%   STYLE = FORMATLATEX(FIG, PRESET) sizes FIG at its intended printed
%   dimensions and applies consistent LaTeX typography to axes, legends,
%   colour bars, tiled-layout titles, and free text.
%
%   Available presets:
%       "half-square"       8.85 x 8.85 cm (0.49 manuscript text width)
%       "full-wide"        18.30 x 6.40 cm (one-by-three layouts)
%       "compact-wide"     15.70 x 7.20 cm (0.86 manuscript text width)
%       "full-landscape"   18.30 x 12.40 cm (two-by-three layouts)
%       "compact-landscape"15.70 x 9.80 cm (0.86 manuscript text width)
%       "full-square"      18.30 x 17.20 cm (two-by-two layouts)
%       "compact-square"   16.50 x 16.50 cm (0.90 manuscript text width)
%
%   Call this function only after all labels, legends, colour bars, and
%   annotations have been created. Plot-specific limits, ticks, grids, and
%   legend locations remain the responsibility of the generating script.

    arguments
        fig (1,1) matlab.ui.Figure
        preset (1,1) string = "full-square"
    end

    style = paperStyle(preset);

    fig.Color = 'w';
    fig.Units = 'centimeters';
    fig.Position(3:4) = [style.Width_cm style.Height_cm];
    fig.PaperUnits = 'centimeters';
    fig.PaperPosition = [0 0 style.Width_cm style.Height_cm];
    fig.PaperSize = [style.Width_cm style.Height_cm];
    fig.PaperPositionMode = 'manual';
    fig.InvertHardcopy = 'off';

    axesHandles = findall(fig, 'Type', 'axes');
    protectedText = gobjects(0);

    for i = 1:numel(axesHandles)
        ax = axesHandles(i);
        ax.TickLabelInterpreter = 'latex';
        ax.FontSize = style.TickFontSize;
        ax.FontWeight = 'normal';
        ax.LineWidth = style.AxesLineWidth;
        ax.Layer = 'top';
        if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
            ax.Toolbar.Visible = 'off';
        end

        labelHandles = [ax.XLabel ax.YLabel ax.ZLabel];
        for j = 1:numel(labelHandles)
            setTextStyle(labelHandles(j), style.LabelFontSize, 'normal');
        end

        % yyaxis stores a separate label on each numeric ruler.
        if isprop(ax, 'YAxis')
            for j = 1:numel(ax.YAxis)
                if isprop(ax.YAxis(j), 'Label')
                    setTextStyle(ax.YAxis(j).Label, style.LabelFontSize, 'normal');
                    labelHandles(end + 1) = ax.YAxis(j).Label; %#ok<AGROW>
                end
            end
        end

        setTextStyle(ax.Title, style.PanelTitleFontSize, 'bold');
        protectedText = [protectedText; labelHandles(:); ax.Title]; %#ok<AGROW>
    end

    legends = findall(fig, 'Type', 'legend');
    for i = 1:numel(legends)
        legends(i).Interpreter = 'latex';
        legends(i).FontSize = style.LegendFontSize;
        legends(i).FontWeight = 'normal';
    end

    colourBars = findall(fig, 'Type', 'colorbar');
    for i = 1:numel(colourBars)
        colourBars(i).TickLabelInterpreter = 'latex';
        colourBars(i).FontSize = style.LegendFontSize;
        if isprop(colourBars(i), 'Label')
            setTextStyle(colourBars(i).Label, style.LabelFontSize, 'normal');
            protectedText(end + 1) = colourBars(i).Label; %#ok<AGROW>
        end
    end

    % Apply the annotation size only to text that is not an axes label or
    % title. Tiled-layout and sgtitle headings are promoted again below.
    freeText = findall(fig, 'Type', 'text');
    if ~isempty(protectedText)
        freeText = setdiff(freeText, protectedText);
    end
    for i = 1:numel(freeText)
        setFontStyle(freeText(i), style.AnnotationFontSize, 'normal');
    end

    textBoxes = findall(fig, 'Type', 'textboxshape');
    for i = 1:numel(textBoxes)
        if isprop(textBoxes(i), 'FontSize')
            textBoxes(i).FontSize = style.AnnotationFontSize;
        end
        if isprop(textBoxes(i), 'FontWeight')
            textBoxes(i).FontWeight = 'normal';
        end
    end

    tiledLayouts = findall(fig, '-isa', 'matlab.graphics.layout.TiledChartLayout');
    for i = 1:numel(tiledLayouts)
        setTextStyle(tiledLayouts(i).Title, style.FigureTitleFontSize, 'bold');
    end

    superTitles = findall(fig, 'Type', 'text', 'Tag', 'suptitle');
    for i = 1:numel(superTitles)
        setTextStyle(superTitles(i), style.FigureTitleFontSize, 'bold');
    end

    drawnow;
end

function style = paperStyle(preset)
    style = struct( ...
        'TickFontSize', 8.0, ...
        'LabelFontSize', 8.5, ...
        'PanelTitleFontSize', 8.5, ...
        'FigureTitleFontSize', 9.5, ...
        'LegendFontSize', 7.25, ...
        'AnnotationFontSize', 7.25, ...
        'AxesLineWidth', 0.65);

    switch lower(preset)
        case "half-square"
            style.Width_cm = 8.85;
            style.Height_cm = 8.85;
        case "full-wide"
            style.Width_cm = 18.30;
            style.Height_cm = 6.40;
        case "compact-wide"
            style.Width_cm = 15.70;
            style.Height_cm = 7.20;
        case "full-landscape"
            style.Width_cm = 18.30;
            style.Height_cm = 12.40;
        case "compact-landscape"
            style.Width_cm = 15.70;
            style.Height_cm = 9.80;
        case "full-square"
            style.Width_cm = 18.30;
            style.Height_cm = 17.20;
        case "compact-square"
            style.Width_cm = 16.50;
            style.Height_cm = 16.50;
        otherwise
            error('formatLatex:UnknownPreset', ...
                'Unknown paper-figure preset "%s".', preset);
    end
end

function setTextStyle(h, fontSize, fontWeight)
    if isempty(h) || ~isgraphics(h)
        return
    end
    if isprop(h, 'Interpreter')
        h.Interpreter = 'latex';
    end
    if isprop(h, 'FontSize')
        h.FontSize = fontSize;
    end
    if isprop(h, 'FontWeight')
        h.FontWeight = fontWeight;
    end
end

function setFontStyle(h, fontSize, fontWeight)
    if isempty(h) || ~isgraphics(h)
        return
    end
    if isprop(h, 'FontSize')
        h.FontSize = fontSize;
    end
    if isprop(h, 'FontWeight')
        h.FontWeight = fontWeight;
    end
end
