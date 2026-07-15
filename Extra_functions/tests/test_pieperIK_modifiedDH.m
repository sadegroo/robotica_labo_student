function tests = test_pieperIK_modifiedDH
% TEST_PIEPERIK_MODIFIEDDH  Unit tests for pieperIK_modifiedDH.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make sure the functions under test are reachable even when the
    % project is not open.  DH_modified is used as the trusted forward
    % map to build target poses and verify IK solutions.
    import matlab.unittest.fixtures.PathFixture
    here = fileparts(mfilename('fullpath'));
    testCase.applyFixture(PathFixture(fullfile(here, '..', 'invkin')));
    testCase.applyFixture(PathFixture(fullfile(here, '..', 'fwdkin')));
    rng(42);  % reproducible random tests
end

function testPUMAExampleRoundTrip(testCase)
    % The PUMA 560 example from the function documentation: the true
    % joint vector must be among the returned solutions.
    DHt = pumaModified();
    theta_true = [10; 20; -30; 40; 50; 60] * pi/180;
    T = fkModified(DHt, theta_true);
    sols = pieperIK_modifiedDH(DHt, T);
    verifyNotEmpty(testCase, sols);
    verifyLessThan(testCase, bestMatch(sols, theta_true), 1e-6, ...
        'true joint vector must be recovered');
end

function testPUMAStatisticalRoundTrip(testCase)
    % 100 random configurations: every solution must reach the pose,
    % the true joint vector must always be recovered, and a generic
    % (non-singular) pose must yield the full set of 8 solutions.
    DHt = pumaModified();
    for trial = 1:100
        theta = (rand(6,1) - 0.5) * 1.4 * pi;
        T = fkModified(DHt, theta);
        sols = pieperIK_modifiedDH(DHt, T);
        verifyEqual(testCase, size(sols), [6 8], ...
            sprintf('trial %d: 8 solutions expected', trial));
        verifyLessThan(testCase, bestMatch(sols, theta), 1e-6, ...
            sprintf('trial %d: true joint vector must be recovered', trial));
        verifyLessThan(testCase, maxPoseError(DHt, sols, T), 1e-9, ...
            sprintf('trial %d: all solutions must reach the pose', trial));
    end
end

function testGenericQuarticBranch(testCase)
    % Robot with a1 ~= 0 and sin(alpha1) ~= 0 exercises the general
    % quartic in tan(theta_3/2).
    DHt = [     0,  0,     0.30,  0;
            -pi/2,  0.15,  0.10,  0;
                0,  0.50,  0.12,  0;
            -pi/2,  0.05,  0.40,  0;
             pi/2,  0,     0,     0;
            -pi/2,  0,     0.08,  0 ];
    roundTripTrials(testCase, DHt, 20);
end

function testNonOrthogonalAlphas(testCase)
    % Twist angles that are not multiples of pi/2, so no term in the
    % derivation degenerates.
    DHt = [     0,  0,     0.30,  0;
             0.50,  0.20,  0.10,  0;
             0.30,  0.45,  0.12,  0;
            -0.40,  0.05,  0.40,  0;
             pi/2,  0,     0,     0;
            -pi/2,  0,     0.08,  0 ];
    roundTripTrials(testCase, DHt, 20);
end

function testAlphaOneZeroBranch(testCase)
    % Robot with sin(alpha1) = 0 and a1 ~= 0 exercises the
    % r_eq-quadratic branch with two (c2, s2) solutions per theta_3.
    DHt = [     0,  0,     0.20,  0;
                0,  0.30,  0.15,  0;
            -pi/2,  0.20,  0.10,  0;
                0,  0.40,  0.30,  0;
            -pi/2,  0,     0,     0;
             pi/2,  0,     0.05,  0 ];
    roundTripTrials(testCase, DHt, 20);
end

function testNonTrivialBaseScrew(testCase)
    % alpha_0 and a_0 nonzero: the solver must pre-transform the target
    % by the inverse base screw and still round-trip.
    DHt = [  0.40,  0.12,  0.30,  0;
            -pi/2,  0.15,  0.10,  0;
             0.30,  0.50,  0.12,  0;
            -pi/2,  0.05,  0.40,  0;
             pi/2,  0,     0,     0;
            -pi/2,  0,     0.08,  0 ];
    roundTripTrials(testCase, DHt, 10);
end

function testSolutionsAreWrapped(testCase)
    % All returned angles must lie in (-pi, pi].
    DHt = pumaModified();
    T = fkModified(DHt, (rand(6,1) - 0.5) * 1.4 * pi);
    sols = pieperIK_modifiedDH(DHt, T);
    verifyNotEmpty(testCase, sols);
    verifyLessThanOrEqual(testCase, max(abs(sols(:))), pi + 1e-12);
end

function testSolutionsAreDistinct(testCase)
    % Away from singularities, the returned solutions must be pairwise
    % distinct under angle wrapping.
    DHt = pumaModified();
    for trial = 1:20
        theta = (rand(6,1) - 0.5) * pi;
        if abs(theta(5)) < 0.2, theta(5) = 0.2 + rand; end
        T = fkModified(DHt, theta);
        sols = pieperIK_modifiedDH(DHt, T);
        for i = 1:size(sols,2)-1
            for j = i+1:size(sols,2)
                d = atan2(sin(sols(:,i) - sols(:,j)), ...
                          cos(sols(:,i) - sols(:,j)));
                verifyGreaterThan(testCase, norm(d), 1e-6, ...
                    sprintf('trial %d: solutions %d and %d duplicate', ...
                            trial, i, j));
            end
        end
    end
end

function testWristSingularityThetaFiveZero(testCase)
    % theta_5 = 0 is a wrist singularity: the solver must warn about it,
    % set theta_4 = 0, and still reach the pose.
    DHt = pumaModified();
    theta = [0.3; -0.4; 0.5; 0.7; 0; 0.9];
    checkSingularPose(testCase, DHt, theta);
end

function testWristSingularityThetaFivePi(testCase)
    % theta_5 = pi is the other wrist singularity (s5 = 0, c5 = -1).
    DHt = pumaModified();
    theta = [0.1; -0.6; 0.4; 1.2; pi; -0.8];
    checkSingularPose(testCase, DHt, theta);
end

function testUnreachablePoseReturnsEmpty(testCase)
    % A target far outside the workspace must give an empty 6x0 result
    % and warn that no real theta_3 roots exist.
    DHt = pumaModified();
    T = [eye(3), [10; 10; 10]; 0 0 0 1];
    lastwarn('');
    evalc('sols = pieperIK_modifiedDH(DHt, T);');    % swallow display
    verifyEqual(testCase, size(sols), [6 0]);
    verifyTrue(testCase, contains(lower(lastwarn()), 'no real'), ...
        'a "no real solutions" warning is expected');
end

function testPieperAssumptionViolationErrors(testCase)
    % Non-zero a_4 = (5,2), a_5 = (6,2) or d_5 = (5,3) must be rejected.
    T = eye(4);
    base = pumaModified();
    violations = {[5 2], [6 2], [5 3]};              % [row col] set to 0.1
    for k = 1:numel(violations)
        DHt = base;
        DHt(violations{k}(1), violations{k}(2)) = 0.1;
        verifyError(testCase, @() pieperIK_modifiedDH(DHt, T), ...
            ?MException, sprintf('violation case %d must error', k));
    end
end

function testInputSizeValidation(testCase)
    verifyError(testCase, ...
        @() pieperIK_modifiedDH(zeros(5,4), eye(4)), ?MException);
    verifyError(testCase, ...
        @() pieperIK_modifiedDH(pumaModified(), eye(3)), ?MException);
end

function testMatchesClassicalSolver(testCase)
    % A classical-DH robot whose trailing x-screw is trivial
    % (alpha_6 = a_6 = 0, as Pieper requires) is EXACTLY the modified-DH
    % robot with the x-screw parameters shifted one link back.  Both
    % solvers must then return the same solution set for the same pose.
    classical = [ -pi/2,  0,       0,       0;
                     0,   0.4318,  0,       0;
                   pi/2,  0.0203,  0.15005, 0;
                  -pi/2,  0,       0.4318,  0;
                   pi/2,  0,       0,       0;
                     0,   0,       0,       0 ];
    modified = classical;
    modified(:,1:2) = [0 0; classical(1:5, 1:2)];    % shift alpha, a
    for trial = 1:10
        theta = (rand(6,1) - 0.5) * pi;
        if abs(theta(5)) < 0.2, theta(5) = 0.2 + rand; end
        T = fk_classicalDH(classical, theta);
        verifyEqual(testCase, fkModified(modified, theta), T, ...
            'AbsTol', 1e-10, 'the two conventions must give the same FK');
        solsC = pieperIK_classicalDH(classical, T);
        solsM = pieperIK_modifiedDH(modified, T);
        verifyEqual(testCase, size(solsM), size(solsC));
        for k = 1:size(solsM, 2)
            verifyLessThan(testCase, bestMatch(solsC, solsM(:,k)), 1e-6, ...
                sprintf('trial %d: modified solution %d missing from classical set', ...
                        trial, k));
        end
    end
end

% =====================================================================
%   Local helpers
% =====================================================================
function DHt = pumaModified()
% Modified-DH PUMA 560 table (Craig ch. 3) from the
% pieperIK_modifiedDH docstring.
    DHt = [    0,   0,       0,       0;
           -pi/2,   0,       0,       0;
               0,   0.4318,  0.15005, 0;
           -pi/2,   0.0203,  0.4318,  0;
            pi/2,   0,       0,       0;
           -pi/2,   0,       0,       0 ];
end

function T = fkModified(DHt, theta)
% Forward kinematics in modified DH via DH_modified.
    tbl = DHt;
    tbl(:,4) = theta;
    T = DH_modified(tbl);
end

function roundTripTrials(testCase, DHt, n_trials)
% Random-configuration round trip: the true joint vector must be
% recovered and every returned solution must reach the target pose.
    for trial = 1:n_trials
        theta = (rand(6,1) - 0.5) * pi;
        if abs(theta(5)) < 0.2, theta(5) = 0.2 + rand; end
        T = fkModified(DHt, theta);
        sols = pieperIK_modifiedDH(DHt, T);
        verifyNotEmpty(testCase, sols, sprintf('trial %d: no solutions', trial));
        verifyLessThanOrEqual(testCase, size(sols, 2), 8, ...
            'at most 8 solutions expected');
        verifyLessThan(testCase, bestMatch(sols, theta), 1e-6, ...
            sprintf('trial %d: true joint vector must be recovered', trial));
        verifyLessThan(testCase, maxPoseError(DHt, sols, T), 1e-9, ...
            sprintf('trial %d: all solutions must reach the pose', trial));
    end
end

function checkSingularPose(testCase, DHt, theta)
% Common body of the wrist-singularity tests: expect a singularity
% warning and solutions that still reach the pose.
    T = fkModified(DHt, theta);
    lastwarn('');
    evalc('sols = pieperIK_modifiedDH(DHt, T);');    % swallow display
    verifyTrue(testCase, contains(lower(lastwarn()), 'singularity'), ...
        'a wrist-singularity warning is expected');
    verifyNotEmpty(testCase, sols);
    verifyLessThan(testCase, maxPoseError(DHt, sols, T), 1e-9);
end

function err = maxPoseError(DHt, sols, T)
% Largest Frobenius-norm pose error over all solution columns.
    err = 0;
    for k = 1:size(sols, 2)
        err = max(err, norm(fkModified(DHt, sols(:,k)) - T, 'fro'));
    end
end

function err = bestMatch(sols, theta_true)
% Smallest max-abs wrapped angular difference between any solution
% column and the true joint vector.
    d = sols - theta_true;
    d = atan2(sin(d), cos(d));          % wrap differences to (-pi, pi]
    err = min(max(abs(d), [], 1));
end
