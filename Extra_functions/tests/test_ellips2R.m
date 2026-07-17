function tests = test_ellips2R
% TEST_ELLIPS2R  Unit tests for ellips2R.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make sure the functions under test are reachable even when the
    % project is not open, and keep test figures off-screen.
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture(fullfile(fileparts(mfilename('fullpath')), '..')));
    rng(42);  % reproducible random tests
    testCase.TestData.oldVisible = get(0, 'DefaultFigureVisible');
    set(0, 'DefaultFigureVisible', 'off');
end

function teardownOnce(testCase)
    set(0, 'DefaultFigureVisible', testCase.TestData.oldVisible);
end

function teardown(testCase) %#ok<INUSD>
    close all
end

function J = jacobian2R(theta, L)
    % Reference spatial 2x2 Jacobian, built independently of ellips2R.
    th1 = theta(1); th12 = theta(1) + theta(2);
    J = [-L(1)*sin(th1) - L(2)*sin(th12), -L(2)*sin(th12);
          L(1)*cos(th1) + L(2)*cos(th12),  L(2)*cos(th12)];
end

function testOutputsMatchManipulabilityEllipsoid(testCase)
    % mel must equal manipulabilityEllipsoid applied to the 2R Jacobian.
    theta = [0.3; 1.1];
    L = [1.2; 0.8];
    [mel, fig] = ellips2R(theta, L);
    expected = manipulabilityEllipsoid(jacobian2R(theta, L));
    verifyEqual(testCase, mel, expected, 'AbsTol', 1e-12);
    verifyClass(testCase, fig, 'matlab.ui.Figure');
end

function testManipulabilityIsUnscaled(testCase)
    % w must be |det(J)| = L1*L2*|sin(th2)|, unaffected by plot scaling.
    theta = [pi/5; 2*pi/3];
    L = [1.5; 0.6];
    mel = ellips2R(theta, L);
    verifyEqual(testCase, mel.w, L(1)*L(2)*abs(sin(theta(2))), 'AbsTol', 1e-12);
end

function testArmPlottedAsTwoSegments(testCase)
    % The arm is one line through base, elbow and end-effector.
    theta = [0.4; 0.9];
    L = [1; 0.7];
    [~, fig] = ellips2R(theta, L);
    ax = findobj(fig, 'Type', 'axes');
    lines = findall(ax, 'Type', 'line');
    arm = lines(arrayfun(@(h) numel(h.XData) == 3, lines));
    verifyNumElements(testCase, arm, 1);
    p_elbow = L(1) * [cos(theta(1)); sin(theta(1))];
    p_ee = p_elbow + L(2) * [cos(sum(theta)); sin(sum(theta))];
    verifyEqual(testCase, [arm.XData; arm.YData], ...
        [[0; 0], p_elbow, p_ee], 'AbsTol', 1e-12);
end

function testEllipseCenteredAtEndEffector(testCase)
    % The mean of the ellipse boundary must be the end-effector position.
    theta = [1.2; -0.7];
    L = [0.9; 1.1];
    [~, fig] = ellips2R(theta, L);
    patchH = findall(fig, 'Tag', 'manipulabilityEllipsoid');
    verifyNumElements(testCase, patchH, 1);
    pts = [patchH.XData(1:end-1).'; patchH.YData(1:end-1).'];  % drop duplicate endpoint
    p_ee = [L(1)*cos(theta(1)) + L(2)*cos(sum(theta));
            L(1)*sin(theta(1)) + L(2)*sin(sum(theta))];
    verifyEqual(testCase, mean(pts, 2), p_ee, 'AbsTol', 1e-10);
end

function testLongestAxisScaledToLinkLengths(testCase)
    % The plotted major semi-axis must be 30% of L1+L2, for any config.
    configs = {[0.3; 1.1], [1.0; 0.2], [-0.5; 2.5]};
    lengths = {[1; 0.7], [2; 0.5], [0.4; 0.4]};
    for k = 1:numel(configs)
        L = lengths{k};
        [~, fig] = ellips2R(configs{k}, L);
        patchH = findall(fig, 'Tag', 'manipulabilityEllipsoid');
        p_ee = [L(1)*cos(configs{k}(1)) + L(2)*cos(sum(configs{k}));
                L(1)*sin(configs{k}(1)) + L(2)*sin(sum(configs{k}))];
        d = hypot(patchH.XData(:) - p_ee(1), patchH.YData(:) - p_ee(2));
        verifyEqual(testCase, max(d), 0.3 * sum(L), 'AbsTol', 1e-10);
    end
end

function testIsotropicConfigurationIsCircle(testCase)
    % L1 = sqrt(2)*L2 and th2 = 3*pi/4 gives a circle with radius L2.
    L = [1; 1/sqrt(2)];
    mel = ellips2R([pi/6; 3*pi/4], L);
    verifyEqual(testCase, mel.sqrt_cd, 1, 'AbsTol', 1e-9);
    verifyEqual(testCase, mel.radii, [L(2); L(2)], 'AbsTol', 1e-9);
end

function testSingularConfigurationDegenerates(testCase)
    % A fully stretched arm has (numerically) zero manipulability.
    mel = ellips2R([pi/6; 0], [1; 0.7]);
    verifyEqual(testCase, mel.w, 0, 'AbsTol', 1e-6);
    verifyEqual(testCase, mel.sqrt_minSV, 0, 'AbsTol', 1e-6);
end

function testTh1OnlyRotatesEllipse(testCase)
    % th1 must not change the manipulability measures, only the pose.
    L = [1.3; 0.9];
    th2 = 0.8;
    melA = ellips2R([0; th2], L);
    melB = ellips2R([2.1; th2], L);
    verifyEqual(testCase, melA.radii, melB.radii, 'AbsTol', 1e-12);
    verifyEqual(testCase, melA.w, melB.w, 'AbsTol', 1e-12);
end

function testInvalidInputsError(testCase)
    verifyError(testCase, @() ellips2R([1, 2], [1; 1]), ?MException);      % row theta
    verifyError(testCase, @() ellips2R([1; 2; 3], [1; 1]), ?MException);   % wrong numel
    verifyError(testCase, @() ellips2R([1; 2], [1, 1]), ?MException);      % row L
    verifyError(testCase, @() ellips2R([1; 2], [1; -1]), ?MException);     % non-positive L
end
