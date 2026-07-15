function tests = test_Euler2RotMat
% TEST_EULER2ROTMAT  Unit tests for Euler2RotMat.
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make sure the functions under test are reachable even when the
    % project is not open.
    import matlab.unittest.fixtures.PathFixture
    here = fileparts(mfilename('fullpath'));
    testCase.applyFixture(PathFixture(fullfile(here, '..', 'fwdkin')));
    rng(42);  % reproducible random tests
end

function testZeroAnglesGiveIdentity(testCase)
    for conv = all_conventions()
        R = Euler2RotMat([0 0 0], conv{1});
        verifyEqual(testCase, R, eye(3), 'AbsTol', 1e-12, ...
            ['identity expected for convention ' conv{1}]);
    end
end

function testKnownZYXEuler(testCase)
    % Compare against explicit composition with the elementary
    % rotation functions of the library.
    a = 0.3; b = -0.5; c = 0.8;
    Rexpected = RotMatZaxis(a, 0) * RotMatYaxis(b, 0) * RotMatXaxis(c, 0);
    R = Euler2RotMat([a b c], 'zyxe');
    verifyEqual(testCase, R, Rexpected, 'AbsTol', 1e-12);
end

function testKnownZYZEuler(testCase)
    % Proper Euler sequence against explicit composition.
    a = 1.1; b = 0.4; c = -2.0;
    Rexpected = RotMatZaxis(a, 0) * RotMatYaxis(b, 0) * RotMatZaxis(c, 0);
    R = Euler2RotMat([a b c], 'zyze');
    verifyEqual(testCase, R, Rexpected, 'AbsTol', 1e-12);
end

function testFixedIsReversedEuler(testCase)
    % 'abcf' must equal the Euler sequence 'cba'e with reversed triple.
    ang = [0.3, -0.7, 1.2];
    for conv = all_conventions()
        c = conv{1};
        if c(4) ~= 'f', continue, end
        Rfixed = Euler2RotMat(ang, c);
        Reuler = Euler2RotMat(ang(end:-1:1), [fliplr(c(1:3)) 'e']);
        verifyEqual(testCase, Rfixed, Reuler, 'AbsTol', 1e-12, ...
            ['fixed/Euler equivalence failed for ' c]);
    end
end

function testFixedComposesInReverseOrder(testCase)
    % 'abcf' gives R = R_c(theta_c) * R_b(theta_b) * R_a(theta_a).
    a = 0.4; b = -1.0; c = 2.1;
    Rexpected = RotMatXaxis(c, 0) * RotMatYaxis(b, 0) * RotMatZaxis(a, 0);
    R = Euler2RotMat([a b c], 'zyxf');
    verifyEqual(testCase, R, Rexpected, 'AbsTol', 1e-12);
end

function testOutputIsInSO3(testCase)
    % Random angles must always produce a proper rotation matrix.
    for conv = all_conventions()
        ang = (rand(1,3) - 0.5) * 2 * pi;
        R = Euler2RotMat(ang, conv{1});
        verifyEqual(testCase, R.'*R, eye(3), 'AbsTol', 1e-10, ...
            'R must be orthogonal');
        verifyEqual(testCase, det(R), 1, 'AbsTol', 1e-10, ...
            'R must have determinant +1');
    end
end

function testConventionIsCaseInsensitive(testCase)
    ang = [0.2, 0.5, -0.9];
    verifyEqual(testCase, Euler2RotMat(ang, 'ZYXE'), ...
        Euler2RotMat(ang, 'zyxe'), 'AbsTol', 1e-14);
end

function testInvalidConventionReturnsNaN(testCase)
    ang = [0.1 0.2 0.3];
    invalid = {'xxze', 'zyxq', 'zyx', 'zyxef', 'abce'};
    for k = 1:numel(invalid)
        R = Euler2RotMat(ang, invalid{k});
        verifyTrue(testCase, all(isnan(R(:))), ...
            ['NaN expected for convention ''' invalid{k} '''']);
    end
end

function testInvalidAnglesReturnNaN(testCase)
    bad = {[1 2], [1 2 3 4], [1 NaN 3], [1 Inf 3], 'abc'};
    for k = 1:numel(bad)
        R = Euler2RotMat(bad{k}, 'zyxe');
        verifyTrue(testCase, all(isnan(R(:))), ...
            sprintf('NaN expected for invalid angles case %d', k));
    end
end

function testSymbolicInput(testCase)
    % Symbolic angles must give a symbolic matrix that matches the
    % numeric result after substitution.
    assumeTrue(testCase, ~isempty(ver('symbolic')), ...
        'Symbolic Math Toolbox not available');
    syms tha thb thc real
    vals = [0.3, -0.7, 1.2];
    for conv = {'zyxe', 'zyze', 'xyzf'}
        Rsym = Euler2RotMat([tha thb thc], conv{1});
        verifyClass(testCase, Rsym, 'sym');
        Rnum = double(subs(Rsym, [tha thb thc], vals));
        verifyEqual(testCase, Rnum, Euler2RotMat(vals, conv{1}), ...
            'AbsTol', 1e-12, ['symbolic mismatch for ' conv{1}]);
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
