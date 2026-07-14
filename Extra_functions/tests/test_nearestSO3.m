function tests = test_nearestSO3
% TEST_NEARESTSO3  Unit tests for nearestSO3.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make sure the function under test is reachable even when the
    % project is not open.
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture(fullfile(fileparts(mfilename('fullpath')), '..')));
    rng(42);  % reproducible random tests
end

function testExactRotationIsFixedPoint(testCase)
    % A matrix already in SO(3) must be returned unchanged with d = 0.
    ax = [1 2 3] / norm([1 2 3]);
    Rin = axang2rotm_local(ax, 0.7);
    [R, d] = nearestSO3(Rin);
    verifyEqual(testCase, R, Rin, 'AbsTol', 1e-12);
    verifyEqual(testCase, d, 0, 'AbsTol', 1e-12);
end

function testIdentity(testCase)
    [R, d] = nearestSO3(eye(3));
    verifyEqual(testCase, R, eye(3), 'AbsTol', 1e-12);
    verifyEqual(testCase, d, 0, 'AbsTol', 1e-12);
end

function testOutputIsInSO3(testCase)
    % For random input the output must be a proper rotation.
    for k = 1:20
        M = randn(3);
        R = nearestSO3(M);
        verifyEqual(testCase, R.'*R, eye(3), 'AbsTol', 1e-10, ...
            'R must be orthogonal');
        verifyEqual(testCase, det(R), 1, 'AbsTol', 1e-10, ...
            'R must have determinant +1');
    end
end

function testDistanceMatchesFrobeniusNorm(testCase)
    M = randn(3);
    [R, d] = nearestSO3(M);
    verifyEqual(testCase, d, norm(M - R, 'fro'), 'AbsTol', 1e-12);
end

function testScaledRotation(testCase)
    % For M = c*Rin (c > 0) the projection is Rin and d = |c-1|*sqrt(3).
    c = 2.5;
    Rin = axang2rotm_local([0 0 1], 1.1);
    [R, d] = nearestSO3(c * Rin);
    verifyEqual(testCase, R, Rin, 'AbsTol', 1e-10);
    verifyEqual(testCase, d, abs(c - 1) * sqrt(3), 'AbsTol', 1e-10);
end

function testNegativeDeterminantInput(testCase)
    % Input with det < 0 (a reflection) must still map to a proper rotation.
    M = diag([1 1 -1]);
    R = nearestSO3(M);
    verifyEqual(testCase, det(R), 1, 'AbsTol', 1e-10);
end

function testOptimality(testCase)
    % No nearby rotation may be closer to M than the returned projection.
    M = randn(3);
    [R, d] = nearestSO3(M);
    for k = 1:50
        ax = randn(1, 3);
        Rother = R * axang2rotm_local(ax / norm(ax), 0.2 * rand);
        verifyLessThanOrEqual(testCase, d, norm(M - Rother, 'fro') + 1e-12);
    end
end

function R = axang2rotm_local(ax, angle)
% Rodrigues formula; avoids a Robotics System Toolbox dependency.
    ax = ax(:) / norm(ax);
    K = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
    R = eye(3) + sin(angle) * K + (1 - cos(angle)) * K^2;
end
