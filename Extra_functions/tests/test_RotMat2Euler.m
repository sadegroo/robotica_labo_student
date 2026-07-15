function tests = test_RotMat2Euler
% TEST_ROTMAT2EULER  Unit tests for RotMat2Euler.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make sure the functions under test are reachable even when the
    % project is not open.  Euler2RotMat is used as the trusted
    % forward map to build test matrices and check round trips.
    import matlab.unittest.fixtures.PathFixture
    here = fileparts(mfilename('fullpath'));
    testCase.applyFixture(PathFixture(fullfile(here, '..', 'invkin')));
    testCase.applyFixture(PathFixture(fullfile(here, '..', 'fwdkin')));
    rng(42);  % reproducible random tests
end

function testKnownZYXAngles(testCase)
    % Angles in the principal range must be recovered exactly as sol1.
    ang = [0.3, -0.5, 0.8];   % beta in (-pi/2, pi/2) -> sol1
    R = Euler2RotMat(ang, 'zyxe');
    [sol1, ~] = RotMat2Euler(R, 'zyxe');
    verifyEqual(testCase, sol1, ang, 'AbsTol', 1e-10);
end

function testKnownZYZAngles(testCase)
    % Proper Euler: beta in (0, pi) is the sol1 branch.
    ang = [1.1, 0.4, -2.0];
    R = Euler2RotMat(ang, 'zyze');
    [sol1, ~] = RotMat2Euler(R, 'zyze');
    verifyEqual(testCase, sol1, ang, 'AbsTol', 1e-10);
end

function testBothSolutionsReproduceR(testCase)
    % Round trip: both returned triples must rebuild the input matrix,
    % for every convention.
    for conv = all_conventions()
        c = conv{1};
        ang = (rand(1,3) - 0.5) * 2 * pi;
        R = Euler2RotMat(ang, c);
        [sol1, sol2] = RotMat2Euler(R, c);
        verifyEqual(testCase, Euler2RotMat(sol1, c), R, 'AbsTol', 1e-9, ...
            ['sol1 round trip failed for ' c]);
        verifyEqual(testCase, Euler2RotMat(sol2, c), R, 'AbsTol', 1e-9, ...
            ['sol2 round trip failed for ' c]);
    end
end

function testSolutionsDifferInGeneral(testCase)
    % Away from gimbal lock the two solutions must be distinct.
    R = Euler2RotMat([0.3, -0.5, 0.8], 'zyxe');
    [sol1, sol2] = RotMat2Euler(R, 'zyxe');
    verifyGreaterThan(testCase, norm(sol1 - sol2), 1e-6);
end

function testGimbalLockTaitBryan(testCase)
    % beta = +/- pi/2 collapses the solutions; theta_a = 0 by convention.
    for beta = [pi/2, -pi/2]
        ang = [0.7, beta, -1.3];
        R = Euler2RotMat(ang, 'zyxe');
        [sol1, sol2] = RotMat2Euler(R, 'zyxe');
        verifyEqual(testCase, sol1, sol2, ...
            'solutions must collapse in gimbal lock');
        verifyEqual(testCase, sol1(1), 0, 'AbsTol', 1e-10, ...
            'theta_a = 0 convention in gimbal lock');
        verifyEqual(testCase, Euler2RotMat(sol1, 'zyxe'), R, ...
            'AbsTol', 1e-9, 'gimbal-lock solution must rebuild R');
    end
end

function testGimbalLockProperEuler(testCase)
    % beta = 0 or pi collapses the solutions for symmetric sequences.
    for beta = [0, pi]
        ang = [0.7, beta, -1.3];
        R = Euler2RotMat(ang, 'zyze');
        [sol1, sol2] = RotMat2Euler(R, 'zyze');
        verifyEqual(testCase, sol1, sol2, ...
            'solutions must collapse in gimbal lock');
        verifyEqual(testCase, sol1(1), 0, 'AbsTol', 1e-10, ...
            'theta_a = 0 convention in gimbal lock');
        verifyEqual(testCase, Euler2RotMat(sol1, 'zyze'), R, ...
            'AbsTol', 1e-9, 'gimbal-lock solution must rebuild R');
    end
end

function testFixedConventionOrdering(testCase)
    % 'abcf' output (a,b,c) must satisfy R = R_c * R_b * R_a, i.e. the
    % reversed triple is the 'cba'e Euler solution.
    ang = [0.3, -0.5, 0.8];
    R = Euler2RotMat(ang, 'xyzf');
    [solF, ~] = RotMat2Euler(R, 'xyzf');
    [solE, ~] = RotMat2Euler(R, 'zyxe');
    verifyEqual(testCase, solF, solE(end:-1:1), 'AbsTol', 1e-10);
end

function testNonRotationInputReturnsNaN(testCase)
    bad = {2*eye(3), diag([1 1 -1]), randn(3)*5, eye(4), zeros(2,3)};
    for k = 1:numel(bad)
        [sol1, sol2] = RotMat2Euler(bad{k}, 'zyxe');
        verifyTrue(testCase, all(isnan([sol1 sol2])), ...
            sprintf('NaN expected for non-rotation input case %d', k));
    end
end

function testInvalidConventionReturnsNaN(testCase)
    invalid = {'xxze', 'zyxq', 'zyx', 'zyxef'};
    for k = 1:numel(invalid)
        [sol1, sol2] = RotMat2Euler(eye(3), invalid{k});
        verifyTrue(testCase, all(isnan([sol1 sol2])), ...
            ['NaN expected for convention ''' invalid{k} '''']);
    end
end

function testSymbolicInput(testCase)
    % A symbolic matrix must yield symbolic extraction formulas that
    % rebuild the matrix after numeric substitution (generic branch).
    assumeTrue(testCase, ~isempty(ver('symbolic')), ...
        'Symbolic Math Toolbox not available');
    syms tha thb thc real
    vals = [0.3, -0.5, 0.8];
    for conv = {'zyxe', 'zyze', 'xyzf'}
        c = conv{1};
        Rsym = Euler2RotMat([tha thb thc], c);
        [sol1sym, ~] = RotMat2Euler(Rsym, c);
        verifyClass(testCase, sol1sym, 'sym');
        angBack = double(subs(sol1sym, [tha thb thc], vals));
        verifyEqual(testCase, Euler2RotMat(angBack, c), ...
            Euler2RotMat(vals, c), 'AbsTol', 1e-10, ...
            ['symbolic extraction failed for ' c]);
    end
end

function convs = all_conventions()
% All 24 valid convention strings (12 sequences x {e,f}).
    convs = {};
    ax = 'xyz';
    for a = 1:3
        for b = 1:3
            for c = 1:3
                if a ~= b && b ~= c
                    for m = 'ef'
                        convs{end+1} = [ax(a) ax(b) ax(c) m]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
