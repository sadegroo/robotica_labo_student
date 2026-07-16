function tests = test_plotEllipsoid
% TEST_PLOTELLIPSOID  Unit tests for plotEllipsoid.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make sure the functions under test are reachable even when the
    % project is not open.
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture(fullfile(fileparts(mfilename('fullpath')), '..')));
    rng(42);  % reproducible random tests
end

function setup(testCase)
    % Every test draws into its own invisible figure.
    testCase.TestData.fig = figure('Visible', 'off');
end

function teardown(testCase)
    close(testCase.TestData.fig);
end

function e = ellipsoid3D()
    e = manipulabilityEllipsoid(randn(3, 5));
end

function e = ellipsoid2D()
    e = manipulabilityEllipsoid(randn(2, 3));
end

function test3DReturnsSurface(testCase)
    h = plotEllipsoid(ellipsoid3D(), [0; 0; 0]);
    verifyClass(testCase, h, 'matlab.graphics.chart.primitive.Surface');
    verifyEqual(testCase, h.Tag, 'manipulabilityEllipsoid');
end

function test2DReturnsPatch(testCase)
    h = plotEllipsoid(ellipsoid2D(), [0; 0]);
    verifyClass(testCase, h, 'matlab.graphics.primitive.Patch');
    verifyEqual(testCase, h.Tag, 'manipulabilityEllipsoid');
end

function test3DSurfaceLiesOnEllipsoid(testCase)
    % Every surface point p must satisfy (p-c)'*inv(A)*(p-c) = 1.
    e = ellipsoid3D();
    c = [0.5; -0.2; 1.0];
    h = plotEllipsoid(e, c);
    P = [h.XData(:).' - c(1); h.YData(:).' - c(2); h.ZData(:).' - c(3)];
    q = sum(P .* (e.A \ P), 1);
    verifyEqual(testCase, q, ones(size(q)), 'AbsTol', 1e-9);
end

function test2DBoundaryLiesOnEllipse(testCase)
    e = ellipsoid2D();
    c = [1; -0.5];
    h = plotEllipsoid(e, c);
    P = [h.XData(:).' - c(1); h.YData(:).' - c(2)];
    q = sum(P .* (e.A \ P), 1);
    verifyEqual(testCase, q, ones(size(q)), 'AbsTol', 1e-9);
end

function testCenterOffsetApplied(testCase)
    % Plotting at center c must translate every point by exactly c
    % compared to the same ellipsoid at the origin.
    e = ellipsoid3D();
    c = [2; 3; -1];
    h0 = plotEllipsoid(e, [0; 0; 0]);
    X0 = h0.XData; Y0 = h0.YData; Z0 = h0.ZData;
    cla(gca, 'reset');
    h1 = plotEllipsoid(e, c);
    verifyEqual(testCase, h1.XData - X0, c(1) * ones(size(X0)), 'AbsTol', 1e-12);
    verifyEqual(testCase, h1.YData - Y0, c(2) * ones(size(Y0)), 'AbsTol', 1e-12);
    verifyEqual(testCase, h1.ZData - Z0, c(3) * ones(size(Z0)), 'AbsTol', 1e-12);
end

function testDefaultCenterIsOrigin(testCase)
    % Omitting the center must equal an explicit zero center.
    e3 = ellipsoid3D();
    hDefault = plotEllipsoid(e3);
    X = hDefault.XData;
    cla(gca, 'reset');
    hZero = plotEllipsoid(e3, [0; 0; 0]);
    verifyEqual(testCase, X, hZero.XData, 'AbsTol', 1e-12);
    cla(gca, 'reset');
    e2 = ellipsoid2D();
    hDefault2 = plotEllipsoid(e2);
    X2 = hDefault2.XData;
    cla(gca, 'reset');
    hZero2 = plotEllipsoid(e2, [0; 0]);
    verifyEqual(testCase, X2, hZero2.XData, 'AbsTol', 1e-12);
end

function testAutoColorsCycle(testCase)
    % Successive calls into the same axes must pick different
    % ColorOrder colors, in both 3D and 2D.
    hold on
    h1 = plotEllipsoid(ellipsoid3D(), [0; 0; 0]);
    h2 = plotEllipsoid(ellipsoid3D(), [1; 0; 0]);
    h3 = plotEllipsoid(ellipsoid3D(), [2; 0; 0]);
    co = colororder(gca);
    verifyEqual(testCase, h1.FaceColor, co(1, :));
    verifyEqual(testCase, h2.FaceColor, co(2, :));
    verifyEqual(testCase, h3.FaceColor, co(3, :));
end

function testAutoColorsCycle2D(testCase)
    hold on
    h1 = plotEllipsoid(ellipsoid2D(), [0; 0]);
    h2 = plotEllipsoid(ellipsoid2D(), [1; 0]);
    co = colororder(gca);
    verifyEqual(testCase, h1.FaceColor, co(1, :));
    verifyEqual(testCase, h2.FaceColor, co(2, :));
end

function testFaceColorOverride(testCase)
    h = plotEllipsoid(ellipsoid3D(), [0; 0; 0], 'FaceColor', [1 0 0]);
    verifyEqual(testCase, h.FaceColor, [1 0 0]);
end

function testFaceAlphaDefaultAndOverride(testCase)
    h = plotEllipsoid(ellipsoid3D(), [0; 0; 0]);
    verifyEqual(testCase, h.FaceAlpha, 0.3);
    cla(gca, 'reset');
    h = plotEllipsoid(ellipsoid3D(), [0; 0; 0], 'FaceAlpha', 0.8);
    verifyEqual(testCase, h.FaceAlpha, 0.8);
end

function testPrincipalAxisLines(testCase)
    % 3 axis lines in 3D, 2 in 2D, colored like their surface and
    % hidden from the legend.
    h = plotEllipsoid(ellipsoid3D(), [0; 0; 0]);
    lines = findall(gca, 'Type', 'line');
    verifyNumElements(testCase, lines, 3);
    for L = lines.'
        verifyEqual(testCase, L.Color, h.FaceColor);
        verifyEqual(testCase, L.HandleVisibility, 'off');
    end
    cla(gca, 'reset');
    h = plotEllipsoid(ellipsoid2D(), [0; 0]);
    lines = findall(gca, 'Type', 'line');
    verifyNumElements(testCase, lines, 2);
    verifyEqual(testCase, lines(1).Color, h.FaceColor);
end

function testAxisLineTipsMatchRadii(testCase)
    % Each principal-axis line must run from the center to
    % center + radii(k)*axes(:,k).
    e = ellipsoid2D();
    c = [1; 2];
    plotEllipsoid(e, c);
    lines = findall(gca, 'Type', 'line');
    pts = [[lines.XData]; [lines.YData]];   % columns: line endpoints
    expected = [c, c + e.radii(1) * e.axes(:, 1), ...
                c, c + e.radii(2) * e.axes(:, 2)];
    verifyEqual(testCase, sortrows(pts.'), sortrows(expected.'), ...
        'AbsTol', 1e-10);
end

function testLegendShowsOneEntryPerEllipsoid(testCase)
    hold on
    plotEllipsoid(ellipsoid3D(), [0; 0; 0], 'DisplayName', 'pose 1');
    plotEllipsoid(ellipsoid3D(), [1; 0; 0], 'DisplayName', 'pose 2');
    lgd = legend('show');
    verifyEqual(testCase, lgd.String, {'pose 1', 'pose 2'});
end

function testHoldStateRestored(testCase)
    % hold off before the call must remain off afterwards, and on
    % must remain on.
    ax = axes(testCase.TestData.fig);
    plotEllipsoid(ellipsoid3D(), [0; 0; 0]);
    verifyFalse(testCase, ishold(ax));
    hold(ax, 'on')
    plotEllipsoid(ellipsoid3D(), [1; 0; 0]);
    verifyTrue(testCase, ishold(ax));
end

function testAxisEqual(testCase)
    plotEllipsoid(ellipsoid3D(), [0; 0; 0]);
    verifyEqual(testCase, get(gca, 'DataAspectRatio'), [1 1 1]);
end

function test2DPlotStaysFlat(testCase)
    % The planar case must produce a genuine 2D plot (default 2D view).
    plotEllipsoid(ellipsoid2D(), [0; 0]);
    [az, el] = view(gca);
    verifyEqual(testCase, [az, el], [0, 90]);
end

function testDegenerateEllipsoidPlots(testCase)
    % A rank-deficient (flat) ellipsoid must plot without error.
    e = manipulabilityEllipsoid([1 2; 2 4]);  % rank 1 -> zero radius
    h = plotEllipsoid(e, [0; 0]);
    verifyClass(testCase, h, 'matlab.graphics.primitive.Patch');
end

function testWrongCenterSizeErrors(testCase)
    verifyError(testCase, ...
        @() plotEllipsoid(ellipsoid3D(), [1; 2]), ?MException);
    verifyError(testCase, ...
        @() plotEllipsoid(ellipsoid2D(), [1; 2; 3]), ?MException);
end
