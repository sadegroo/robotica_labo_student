function h = plotEllipsoid(ellipsoid, center, varargin)
% PLOTELLIPSOID  Plot an ellipsoid from manipulabilityEllipsoid data.
%   h = plotEllipsoid(ellipsoid, center) draws the ellipsoid described by
%   the struct returned by manipulabilityEllipsoid (fields 'axes' and
%   'radii'), centered at center, into the current axes. For 3D data
%   (3 radii) it draws a translucent surface and center must be a
%   3-vector; for planar data (2 radii) it draws a flat filled ellipse
%   in a regular 2D plot and center must be a 2-vector. Returns the
%   surface/patch handle h.
%
%   Each call picks the next color of the axes ColorOrder, so overlaid
%   ellipsoids in the same axes are automatically distinguishable. The
%   principal-axis lines are drawn in the same color as the surface.
%
%   h = plotEllipsoid(ellipsoid, center, Name, Value, ...) forwards extra
%   name/value pairs to surf (3D) or patch (2D). Useful ones:
%     'FaceColor'   - override the automatic color, e.g. 'r' or [0 .5 1]
%     'FaceAlpha'   - transparency, default 0.3
%     'DisplayName' - label shown by legend (one entry per ellipsoid)
%
%   Example (two configurations overlaid):
%     hold on
%     plotEllipsoid(e1, c1, 'DisplayName', 'pose 1');
%     plotEllipsoid(e2, c2, 'DisplayName', 'pose 2');
%     legend show
%
%   See also manipulabilityEllipsoid.

    dim = numel(ellipsoid.radii);
    if nargin < 2 || isempty(center)
        center = zeros(dim, 1);
    end
    validateattributes(center, {'numeric'}, {'vector', 'numel', dim});
    center = double(center(:));

    ax = newplot;

    % Pick the next ColorOrder color unless the caller set FaceColor
    names = varargin(1:2:end);
    idx = find(strcmpi(names, 'FaceColor'), 1, 'last');
    if isempty(idx)
        nPrev = numel(findobj(ax, 'Tag', 'manipulabilityEllipsoid'));
        co = colororder(ax);
        faceColor = co(mod(nPrev, size(co, 1)) + 1, :);
        varargin = [{'FaceColor', faceColor}, varargin];
    else
        faceColor = varargin{2 * idx};
    end

    if dim == 2
        % Planar case: unit circle mapped through the axes and radii
        t = linspace(0, 2*pi, 100);
        pts = ellipsoid.axes * diag(ellipsoid.radii) * [cos(t); sin(t)];
        h = patch(ax, 'XData', pts(1, :) + center(1), ...
                  'YData', pts(2, :) + center(2), ...
                  'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
                  'Tag', 'manipulabilityEllipsoid', varargin{:});
    else
        % Unit sphere mapped through the principal axes and radii
        nGrid = 40;
        [xs, ys, zs] = sphere(nGrid);
        pts = ellipsoid.axes * diag(ellipsoid.radii) * ...
              [xs(:), ys(:), zs(:)].';
        X = reshape(pts(1, :), size(xs)) + center(1);
        Y = reshape(pts(2, :), size(ys)) + center(2);
        Z = reshape(pts(3, :), size(zs)) + center(3);
        h = surf(ax, X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
                 'Tag', 'manipulabilityEllipsoid', varargin{:});
    end

    % Principal semi-axes in the surface color; hidden from the legend
    if ischar(faceColor) || isstring(faceColor)
        if any(strcmpi(faceColor, {'none', 'flat', 'interp', 'texturemap'}))
            lineColor = 'k';
        else
            lineColor = faceColor;
        end
    else
        lineColor = faceColor;
    end
    washold = ishold(ax);
    hold(ax, 'on')
    for k = 1:dim
        a = ellipsoid.radii(k) * ellipsoid.axes(:, k);
        if dim == 2
            plot(ax, center(1) + [0, a(1)], center(2) + [0, a(2)], ...
                 '-', 'Color', lineColor, 'LineWidth', 1, ...
                 'HandleVisibility', 'off');
        else
            plot3(ax, center(1) + [0, a(1)], center(2) + [0, a(2)], ...
                  center(3) + [0, a(3)], '-', 'Color', lineColor, ...
                  'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
    if ~washold
        hold(ax, 'off')
    end
    axis(ax, 'equal')
end
