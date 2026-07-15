function tests = test_DH_modified
% TEST_DH_MODIFIED  Unit tests for DH_modified (Craig's modified DH).
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

function testZeroRowGivesIdentity(testCase)
    T = DH_modified([0 0 0 0]);
    verifyEqual(testCase, T, eye(4), 'AbsTol', 1e-12);
end

function testSingleRowMatchesCraigFormula(testCase)
    % Compare against Craig eq. 3.6 written out explicitly.
    alpha = 0.7; a = 0.3; d = -0.2; theta = 1.1;
    ct = cos(theta); st = sin(theta);
    ca = cos(alpha); sa = sin(alpha);
    Texpected = [ ct,    -st,     0,    a;
                  st*ca,  ct*ca, -sa,  -d*sa;
                  st*sa,  ct*sa,  ca,   d*ca;
                  0,      0,      0,    1 ];
    T = DH_modified([alpha a d theta]);
    verifyEqual(testCase, T, Texpected, 'AbsTol', 1e-12);
end

function testSingleRowIsScrewComposition(testCase)
    % ^{i-1}T_i must equal Rot(x,alpha)*Trans(x,a)*Rot(z,theta)*Trans(z,d).
    alpha = -1.3; a = 0.25; d = 0.4; theta = 0.9;
    Tx = eye(4); Tx(1,4) = a;
    Tz = eye(4); Tz(3,4) = d;
    Texpected = RotMatXaxis(alpha, 1) * Tx * RotMatZaxis(theta, 1) * Tz;
    T = DH_modified([alpha a d theta]);
    verifyEqual(testCase, T, Texpected, 'AbsTol', 1e-12);
end

function testPureThetaIsZRotation(testCase)
    theta = 0.6;
    verifyEqual(testCase, DH_modified([0 0 0 theta]), ...
        RotMatZaxis(theta, 1), 'AbsTol', 1e-12);
end

function testPureAlphaIsXRotation(testCase)
    alpha = -0.8;
    verifyEqual(testCase, DH_modified([alpha 0 0 0]), ...
        RotMatXaxis(alpha, 1), 'AbsTol', 1e-12);
end

function testPureTranslations(testCase)
    a = 0.5; d = 0.3;
    T = DH_modified([0 a d 0]);
    verifyEqual(testCase, T(1:3,1:3), eye(3), 'AbsTol', 1e-12);
    verifyEqual(testCase, T(1:3,4), [a; 0; d], 'AbsTol', 1e-12);
end

function testPartsAndCumulConsistency(testCase)
    % Tcumul must be the running product of Tparts and end at Tfull.
    table = randomTable(5);
    [Tfull, Tparts, Tcumul] = DH_modified(table);
    verifyEqual(testCase, numel(Tparts), 5);
    verifyEqual(testCase, numel(Tcumul), 5);
    running = eye(4);
    for i = 1:5
        running = running * Tparts{i};
        verifyEqual(testCase, Tcumul{i}, running, 'AbsTol', 1e-12, ...
            sprintf('Tcumul{%d} mismatch', i));
    end
    verifyEqual(testCase, Tfull, Tcumul{end}, 'AbsTol', 1e-12);
end

function testRotationPartIsInSO3(testCase)
    % The rotation block of every transform must be a proper rotation.
    [~, Tparts, Tcumul] = DH_modified(randomTable(6));
    for i = 1:6
        for T = {Tparts{i}, Tcumul{i}}
            R = T{1}(1:3,1:3);
            verifyEqual(testCase, R.'*R, eye(3), 'AbsTol', 1e-10);
            verifyEqual(testCase, det(R), 1, 'AbsTol', 1e-10);
            verifyEqual(testCase, T{1}(4,:), [0 0 0 1], 'AbsTol', 1e-12);
        end
    end
end

function testEquivalenceWithClassicalDH(testCase)
    % A chain in classical DH (alpha_i, a_i, d_i, theta_i) equals the
    % modified-DH chain with the x-screw parameters shifted one link
    % back (alpha_0 = a_0 = 0), up to the trailing x-screw of link n:
    %   T_classical = T_modified * Rot(x, alpha_n) * Trans(x, a_n)
    n = 6;
    classical = randomTable(n);
    modified  = classical;
    modified(:,1:2)     = [0 0; classical(1:n-1, 1:2)];  % shift alpha, a
    Tclass = DH_full(classical);
    Tmod   = DH_modified(modified);
    Tx = eye(4); Tx(1,4) = classical(n,2);
    verifyEqual(testCase, Tmod * RotMatXaxis(classical(n,1), 1) * Tx, ...
        Tclass, 'AbsTol', 1e-10);
end

function testSymbolicInput(testCase)
    % A symbolic table must give symbolic transforms matching the
    % numeric result after substitution.
    assumeTrue(testCase, ~isempty(ver('symbolic')), ...
        'Symbolic Math Toolbox not available');
    syms al a d th real
    [Tfull, Tparts, Tcumul] = DH_modified([al a d th; -al 2*a -d th/2]);
    verifyClass(testCase, Tfull, 'sym');
    vals = [0.7, 0.3, -0.2, 1.1];
    table_num = [vals; -vals(1), 2*vals(2), -vals(3), vals(4)/2];
    [Tfull_num, Tparts_num, Tcumul_num] = DH_modified(table_num);
    verifyEqual(testCase, double(subs(Tfull, [al a d th], vals)), ...
        Tfull_num, 'AbsTol', 1e-12);
    for i = 1:2
        verifyEqual(testCase, double(subs(Tparts{i}, [al a d th], vals)), ...
            Tparts_num{i}, 'AbsTol', 1e-12);
        verifyEqual(testCase, double(subs(Tcumul{i}, [al a d th], vals)), ...
            Tcumul_num{i}, 'AbsTol', 1e-12);
    end
end

% =====================================================================
%   Local helpers
% =====================================================================
function table = randomTable(n)
% Random n x 4 modified/classical DH table with angles in (-pi, pi)
% and link lengths/offsets in (-0.5, 0.5).
    table = [(rand(n,1) - 0.5) * 2*pi, ...
             (rand(n,2) - 0.5), ...
             (rand(n,1) - 0.5) * 2*pi];
end
