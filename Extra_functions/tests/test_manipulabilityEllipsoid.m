function tests = test_manipulabilityEllipsoid
% TEST_MANIPULABILITYELLIPSOID  Unit tests for manipulabilityEllipsoid.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make sure the function under test is reachable even when the
    % project is not open.
    import matlab.unittest.fixtures.PathFixture
    testCase.applyFixture(PathFixture(fullfile(fileparts(mfilename('fullpath')), '..')));
    rng(42);  % reproducible random tests
end

function testLinearIsDefaultFor6xn(testCase)
    % Default on a 6xn Jacobian must equal the explicit 'linear' choice
    % and use the bottom three rows.
    J = randn(6, 4);
    eDefault = manipulabilityEllipsoid(J);
    eLinear = manipulabilityEllipsoid(J, 'linear');
    verifyEqual(testCase, eDefault, eLinear);
    verifyEqual(testCase, eDefault.A, J(4:6, :) * J(4:6, :).', 'AbsTol', 1e-12);
end

function testAngularUsesTopRows(testCase)
    J = randn(6, 5);
    e = manipulabilityEllipsoid(J, 'angular');
    verifyEqual(testCase, e.A, J(1:3, :) * J(1:3, :).', 'AbsTol', 1e-12);
end

function testPartIsCaseInsensitive(testCase)
    J = randn(6, 3);
    verifyEqual(testCase, manipulabilityEllipsoid(J, 'Angular'), ...
        manipulabilityEllipsoid(J, 'angular'));
end

function test3xnUsedAsIs(testCase)
    % A 3xn Jacobian is used directly and 'part' is ignored.
    J = randn(3, 4);
    e = manipulabilityEllipsoid(J);
    verifyEqual(testCase, e.A, J * J.', 'AbsTol', 1e-12);
    verifyEqual(testCase, manipulabilityEllipsoid(J, 'angular'), e);
end

function test2xnPlanarGives2DEllipse(testCase)
    % A planar Jacobian must yield 2x2 axes, 2x1 radii, 2x2 A.
    J = randn(2, 3);
    e = manipulabilityEllipsoid(J);
    verifySize(testCase, e.axes, [2 2]);
    verifySize(testCase, e.radii, [2 1]);
    verifySize(testCase, e.A, [2 2]);
    verifyEqual(testCase, e.A, J * J.', 'AbsTol', 1e-12);
end

function testRadiiAreSingularValues(testCase)
    % Semi-axis lengths equal the singular values of the used block,
    % sorted in descending order, for every input shape.
    Js = {randn(6, 4), randn(3, 5), randn(2, 2)};
    blocks = {@(J) J(4:6, :), @(J) J, @(J) J};
    for k = 1:numel(Js)
        e = manipulabilityEllipsoid(Js{k});
        sv = sort(svd(blocks{k}(Js{k})), 'descend');
        verifyEqual(testCase, e.radii, sv, 'AbsTol', 1e-10);
        verifyTrue(testCase, issorted(e.radii, 'descend'));
    end
end

function testAxesAreRightHandedOrthonormal(testCase)
    for J = {randn(6, 6), randn(3, 4), randn(2, 3)}
        e = manipulabilityEllipsoid(J{1});
        n = size(e.axes, 1);
        verifyEqual(testCase, e.axes.' * e.axes, eye(n), 'AbsTol', 1e-10, ...
            'principal axes must be orthonormal');
        verifyEqual(testCase, det(e.axes), 1, 'AbsTol', 1e-10, ...
            'principal axes must form a right-handed set');
    end
end

function testEigendecompositionReconstructsA(testCase)
    % axes * diag(radii^2) * axes' must reproduce A.
    J = randn(6, 5);
    e = manipulabilityEllipsoid(J);
    verifyEqual(testCase, e.axes * diag(e.radii.^2) * e.axes.', e.A, ...
        'AbsTol', 1e-10);
end

function testManipulabilityMeasure(testCase)
    % w = prod(radii) = sqrt(det(A)).
    J = randn(3, 6);
    e = manipulabilityEllipsoid(J);
    verifyEqual(testCase, e.w, sqrt(det(J * J.')), 'AbsTol', 1e-10);
    verifyEqual(testCase, e.w, prod(e.radii), 'AbsTol', 1e-12);
end

function testJointVelocityImageInsideEllipsoid(testCase)
    % The ellipsoid is the image of the unit sphere in joint space:
    % every Jp*qd with ||qd|| = 1 must satisfy p'*inv(A)*p <= 1, and
    % the principal directions are reached exactly.
    J = randn(3, 5);
    e = manipulabilityEllipsoid(J);
    for k = 1:25
        qd = randn(5, 1); qd = qd / norm(qd);
        p = J * qd;
        verifyLessThanOrEqual(testCase, p.' * (e.A \ p), 1 + 1e-9);
    end
end

function testRankDeficientJacobian(testCase)
    % A singular configuration must give a zero smallest radius and
    % real, non-negative output (no complex numbers from round-off).
    J = [1 2; 2 4];  % rank 1
    e = manipulabilityEllipsoid(J);
    verifyTrue(testCase, isreal(e.radii) && isreal(e.axes));
    verifyGreaterThanOrEqual(testCase, e.radii, 0);
    verifyEqual(testCase, e.radii(2), 0, 'AbsTol', 1e-10);
    verifyEqual(testCase, e.w, 0, 'AbsTol', 1e-10);
end

function testZeroJacobian(testCase)
    e = manipulabilityEllipsoid(zeros(3, 4));
    verifyEqual(testCase, e.radii, zeros(3, 1), 'AbsTol', 1e-12);
    verifyEqual(testCase, e.w, 0);
end

function testKnownDiagonalCase(testCase)
    % For Jp = diag([3 2 1]) the ellipsoid axes align with x, y, z and
    % the radii are exactly 3, 2, 1.
    e = manipulabilityEllipsoid([zeros(3); diag([3 2 1])]);
    verifyEqual(testCase, e.radii, [3; 2; 1], 'AbsTol', 1e-12);
    verifyEqual(testCase, abs(e.axes), eye(3), 'AbsTol', 1e-12);
end

function testInvalidPartErrors(testCase)
    verifyError(testCase, @() manipulabilityEllipsoid(randn(6, 3), 'lineair'), ...
        ?MException);
end

function testInvalidRowCountErrors(testCase)
    verifyError(testCase, @() manipulabilityEllipsoid(randn(4, 4)), ...
        ?MException);
    verifyError(testCase, @() manipulabilityEllipsoid(randn(1, 4)), ...
        ?MException);
end
