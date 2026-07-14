function tests = test_nearestSE3
% TEST_NEARESTSE3  Unit tests for nearestSE3.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make sure the functions under test are reachable even when the
    % project is not open.
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture(fullfile(fileparts(mfilename('fullpath')), '..')));
    rng(42);  % reproducible random tests
end

function testExactSE3IsFixedPoint(testCase)
    % A matrix already in SE(3) must be returned unchanged with d = 0.
    Tin = makeSE3([0.4 -0.8 1.2], [1 2 3]);
    [T, d] = nearestSE3(Tin);
    verifyEqual(testCase, T, Tin, 'AbsTol', 1e-12);
    verifyEqual(testCase, d, 0, 'AbsTol', 1e-12);
end

function testOutputStructure(testCase)
    % For random input: rotation block in SO(3), translation preserved,
    % bottom row exactly [0 0 0 1].
    for k = 1:20
        Ttilde = randn(4);
        T = nearestSE3(Ttilde);
        R = T(1:3, 1:3);
        verifyEqual(testCase, R.'*R, eye(3), 'AbsTol', 1e-10);
        verifyEqual(testCase, det(R), 1, 'AbsTol', 1e-10);
        verifyEqual(testCase, T(1:3, 4), Ttilde(1:3, 4), ...
            'Translation must be copied unchanged');
        verifyEqual(testCase, T(4, :), [0 0 0 1], ...
            'Bottom row must be exactly [0 0 0 1]');
    end
end

function testDistanceMatchesFrobeniusNorm(testCase)
    % Because the translation is preserved exactly, d must equal the full
    % Frobenius distance ||Ttilde - T||_F.
    Ttilde = randn(4);
    [T, d] = nearestSE3(Ttilde);
    verifyEqual(testCase, d, norm(Ttilde - T, 'fro'), 'AbsTol', 1e-12);
end

function testNoisyTransformRecovered(testCase)
    % Small perturbation of a valid transform must project back close to it.
    Tin = makeSE3([0.2 0.5 -0.3], [-1 0.5 2]);
    noise = 1e-6 * randn(4);
    [T, d] = nearestSE3(Tin + noise);
    verifyEqual(testCase, T, Tin, 'AbsTol', 1e-4);
    verifyLessThan(testCase, d, 1e-4);
end

function testRotationBlockMatchesNearestSO3(testCase)
    % The rotation block must be exactly the nearestSO3 projection of the
    % upper-left 3x3 block.
    Ttilde = randn(4);
    T = nearestSE3(Ttilde);
    Rexpected = nearestSO3(Ttilde(1:3, 1:3));
    verifyEqual(testCase, T(1:3, 1:3), Rexpected, 'AbsTol', 1e-12);
end

function testBottomRowPenalty(testCase)
    % A drifted bottom row contributes exactly its distance to [0 0 0 1].
    Tin = makeSE3([0 0 0.9], [1 -2 3]);
    Ttilde = Tin;
    Ttilde(4, :) = [0.1 -0.2 0.3 1.4];
    [T, d] = nearestSE3(Ttilde);
    verifyEqual(testCase, T, Tin, 'AbsTol', 1e-12);
    verifyEqual(testCase, d, norm([0.1 -0.2 0.3 0.4]), 'AbsTol', 1e-12);
end

function T = makeSE3(rotvec, t)
% Build a valid SE(3) matrix from a rotation vector and translation.
    angle = norm(rotvec);
    if angle == 0
        R = eye(3);
    else
        ax = rotvec(:) / angle;
        K = [0 -ax(3) ax(2); ax(3) 0 -ax(1); -ax(2) ax(1) 0];
        R = eye(3) + sin(angle) * K + (1 - cos(angle)) * K^2;
    end
    T = [R, t(:); 0 0 0 1];
end
