function reserveLegendMetricMargin(lgd, marginFraction, boundsAxes)
%RESERVELEGENDMETRICMARGIN Reserve width for LaTeX text in vector exports.
%   RESERVELEGENDMETRICMARGIN(LGD) widens a boxed legend by 8 percent after
%   its final font and placement have been resolved. This compensates for
%   the slightly wider glyph metrics used by MATLAB's vector PDF renderer.
%
%   RESERVELEGENDMETRICMARGIN(LGD, FRACTION) uses a custom fractional width
%   allowance. Call this function after FORMATLATEX and before exporting.
%
%   RESERVELEGENDMETRICMARGIN(LGD, FRACTION, AX) additionally constrains the
%   widened legend to the plotting rectangle of axes AX. This is useful for
%   legends in tiled layouts, where automatic placement can extend into the
%   gap between neighbouring panels.

    if nargin < 2
        marginFraction = 0.08;
    end
    if nargin < 3
        boundsAxes = gobjects(0);
    end
    validateattributes(marginFraction, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'}, ...
        mfilename, 'marginFraction');

    if isempty(lgd) || ~isgraphics(lgd) || strcmpi(lgd.Box, 'off')
        return
    end

    drawnow;
    originalUnits = lgd.Units;
    originalLocation = lower(string(lgd.Location));
    lgd.Units = 'normalized';
    pos = lgd.Position;

    extraWidth = marginFraction * pos(3);
    if contains(originalLocation, "east")
        preserveRightEdge = true;
    elseif contains(originalLocation, "west")
        preserveRightEdge = false;
    else
        % For automatic or centre-based locations, infer the anchor from
        % the resolved position rather than the unresolved location name.
        preserveRightEdge = pos(1) + 0.5 * pos(3) >= 0.5;
    end

    if preserveRightEdge
        pos(1) = pos(1) - extraWidth;
    end
    pos(3) = pos(3) + extraWidth;

    if ~isempty(boundsAxes) && isgraphics(boundsAxes, 'axes')
        originalAxesUnits = boundsAxes.Units;
        boundsAxes.Units = 'normalized';
        axesPos = boundsAxes.Position;
        boundsAxes.Units = originalAxesUnits;

        horizontalInset = 0.01 * axesPos(3);
        minimumX = axesPos(1) + horizontalInset;
        maximumX = axesPos(1) + axesPos(3) - horizontalInset;
        availableWidth = maximumX - minimumX;
        pos(3) = min(pos(3), availableWidth);
        pos(1) = max(minimumX, min(pos(1), maximumX - pos(3)));
    else
        % Fall back to the fixed figure canvas when no axes bound is given.
        pos(1) = max(0.01, min(pos(1), 0.99 - pos(3)));
    end
    lgd.Position = pos;
    lgd.AutoUpdate = 'off';
    lgd.Units = originalUnits;
    drawnow;
end
