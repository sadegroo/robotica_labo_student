function [mel, fig] = ellips2R(theta, L)
% ELLIPS2R  Plot de manipuleerbaarheidsellips van een planaire 2R arm.
%   [mel, fig] = ellips2R(theta, L) plot de arm (2 lijnsegmenten) en de
%   manipuleerbaarheidsellips, gecentreerd op de eindeffector. De ellips
%   wordt geschaald zodat de langste halve as steeds 30% van L1+L2
%   bedraagt.
%
%   Inputs:
%     theta - 2x1 kolomvector gewrichtshoeken [th1; th2] (rad)
%     L     - 2x1 kolomvector schakel-lengtes [L1; L2]
%   Outputs:
%     mel - struct van manipulabilityEllipsoid met de (ongeschaalde)
%           manipuleerbaarheidswaarden: w, sqrt_minSV, sqrt_cd, ...
%     fig - handle van de figuur

    validateattributes(theta, {'numeric'}, {'column', 'numel', 2});
    validateattributes(L, {'numeric'}, {'column', 'numel', 2, 'positive'});

    th1 = theta(1);
    th12 = theta(1) + theta(2);
    L1 = L(1);
    L2 = L(2);

    % Voorwaartse kinematica: positie van elleboog en eindeffector
    p_elbow = [L1*cos(th1); L1*sin(th1)];
    p_ee = [L1*cos(th1) + L2*cos(th12); L1*sin(th1) + L2*sin(th12)];

    % Ruimtelijke 2x2 Jacobiaan: afgeleide van de kinematica naar theta
    J = [-L1*sin(th1) - L2*sin(th12), -L2*sin(th12);
          L1*cos(th1) + L2*cos(th12),  L2*cos(th12)];

    mel = manipulabilityEllipsoid(J);

    % Schaal voor het plotten: langste halve as = 30% van L1+L2
    melPlot = mel;
    melPlot.radii = mel.radii * (0.3 * (L1 + L2) / mel.radii(1));

    fig = figure;
    plot([0, p_elbow(1), p_ee(1)], [0, p_elbow(2), p_ee(2)], '-o', ...
         'LineWidth', 2, 'Color', [0.2 0.2 0.2], 'MarkerFaceColor', 'w');
    hold on
    plotEllipsoid(melPlot, p_ee);
    hold off
    grid on
    % Krappe aslimieten: bounding box van basis, elleboog en ellips
    ext = sqrt(sum((melPlot.axes * diag(melPlot.radii)).^2, 2));
    pts = [[0; 0], p_elbow, p_ee - ext, p_ee + ext];
    pad = 0.1 * (L1 + L2);
    xlim([min(pts(1, :)) - pad, max(pts(1, :)) + pad]);
    ylim([min(pts(2, :)) - pad, max(pts(2, :)) + pad]);
    xlabel('x');
    ylabel('y');
    title(sprintf('2R arm: \\theta = [%.2f; %.2f] rad', theta(1), theta(2)));
    subtitle(sprintf('w = %.3g,  sqrt_minSV = %.3g,  sqrt_cd = %.3g', ...
                     mel.w, mel.sqrt_minSV, mel.sqrt_cd), ...
             'Interpreter', 'none');
end