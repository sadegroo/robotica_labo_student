function [sol1, sol2] = RotMat2Euler(R, convention)
%RotMat2Euler  Extract Euler or fixed-angle triples from a rotation matrix.
%
%   [sol1, sol2] = RotMat2Angles(R, convention) returns both
%   solutions for the chosen convention.  Each solution is a 1x3 row
%   vector of angles in radians.
%
%   convention is a 4-character string:
%       characters 1-3 : axis sequence, each in {'x','y','z'}, with
%                        adjacent axes different (12 valid sequences)
%       character 4    : 'e' for Euler (rotations about the moving frame)
%                        'f' for fixed (rotations about the fixed frame)
%   Combined, this covers all 24 conventions of Craig's Appendix B.
%
%   Output ordering:
%       'abce' returns [theta_a, theta_b, theta_c] such that
%               R = R_a(theta_a) * R_b(theta_b) * R_c(theta_c).
%       'abcf' returns [theta_a, theta_b, theta_c] such that
%               R = R_c(theta_c) * R_b(theta_b) * R_a(theta_a).
%
%   In gimbal lock (cos(beta) = 0 for Tait-Bryan, sin(beta) = 0 for
%   proper Euler) the two solutions collapse and the convention
%   theta_a = 0 is applied; sol2 is returned equal to sol1.
%
%   If R is not a proper rotation matrix or the convention string is
%   invalid, both outputs are returned as NaN(1,3).
%
%   R may be symbolic (sym).  In that case orthonormality is not
%   checked and the generic (non-gimbal-lock) formulas are returned,
%   since the degenerate branch cannot be decided symbolically.

    % --- Validate R -----------------------------------------------------
    tol = 1e-6;
    if ~isequal(size(R), [3 3])
        sol1 = nan(1,3); sol2 = nan(1,3); return
    end
    if ~isa(R, 'sym') && ...
            (norm(R.'*R - eye(3), 'fro') > tol || abs(det(R) - 1) > tol)
        sol1 = nan(1,3); sol2 = nan(1,3); return
    end

    % --- Validate convention --------------------------------------------
    convention = lower(char(convention));
    if numel(convention) ~= 4 || ...
            ~all(ismember(convention(1:3), 'xyz')) || ...
            ~ismember(convention(4), {'e','f'}) || ...
            convention(1) == convention(2) || ...
            convention(2) == convention(3)
        sol1 = nan(1,3); sol2 = nan(1,3); return
    end

    axes_str = convention(1:3);
    mode     = convention(4);

    % --- Map fixed -> equivalent Euler ----------------------------------
    %   'abcf' with output (g,b,a) == 'cba'e with output (a,b,g),
    %   so we extract Euler with the reversed sequence and then
    %   reverse the output triple.
    if mode == 'f'
        euler_axes = fliplr(axes_str);
    else
        euler_axes = axes_str;
    end

    ax = struct('x',1, 'y',2, 'z',3);
    i = ax.(euler_axes(1));
    j = ax.(euler_axes(2));
    k = ax.(euler_axes(3));

    if i == k
        % Proper Euler (symmetric, e.g. Z-Y-Z)
        l     = setdiff([1 2 3], [i j]);    % the unused axis
        sigma = perm_sign([i j l]);
        [s1, s2] = solve_proper(R, i, j, l, sigma, tol);
    else
        % Tait-Bryan (asymmetric, e.g. Z-Y-X)
        sigma = perm_sign([i j k]);
        [s1, s2] = solve_taitbryan(R, i, j, k, sigma, tol);
    end

    if mode == 'f'
        sol1 = fliplr(s1);
        sol2 = fliplr(s2);
    else
        sol1 = s1;
        sol2 = s2;
    end
end


% =====================================================================
%   Tait-Bryan extraction:  R = R_i(a) * R_j(b) * R_k(c)
%   Generic identities (sigma = +/- 1 from the permutation parity):
%       sin(b)            =  sigma * R(i,k)
%       cos(b)            =  +/- sqrt( R(i,i)^2 + R(i,j)^2 )
%       tan(a)            =  -sigma * R(j,k) / R(k,k)
%       tan(c)            =  -sigma * R(i,j) / R(i,i)
% =====================================================================
function [s1, s2] = solve_taitbryan(R, i, j, k, sigma, tol)
    sin_b = sigma * R(i,k);
    if ~isa(R, 'sym')
        sin_b = max(min(sin_b, 1), -1);         % clamp against round-off
    end
    cb_sq = R(i,i)^2 + R(i,j)^2;

    if isa(R, 'sym') || cb_sq > tol^2
        cos_b = sqrt(cb_sq);
        b1 = atan2(sin_b,  cos_b);              % b1 in (-pi/2, pi/2)
        a1 = atan2(-sigma*R(j,k),  R(k,k));
        c1 = atan2(-sigma*R(i,j),  R(i,i));

        b2 = atan2(sin_b, -cos_b);              % b2 = pi - b1  (mod 2pi)
        a2 = atan2( sigma*R(j,k), -R(k,k));
        c2 = atan2( sigma*R(i,j), -R(i,i));

        s1 = [a1, b1, c1];
        s2 = [a2, b2, c2];
    else
        % Gimbal lock: b = +/- pi/2.  Set a = 0 and read c off the
        % residual matrix R_j(-b) * R, which collapses to R_k(c).
        b1 = atan2(sin_b, 0);                   % +pi/2 or -pi/2
        a1 = 0;
        residual = elem_rot(j, -b1) * R;
        c1 = extract_elementary(residual, k);
        s1 = [a1, b1, c1];
        s2 = s1;
    end
end


% =====================================================================
%   Proper Euler extraction:  R = R_i(a) * R_j(b) * R_i(c)
%   Letting l be the unused axis and sigma the parity of (i,j,l):
%       cos(b)            =  R(i,i)
%       sin(b)            =  +/- sqrt( R(i,j)^2 + R(i,l)^2 )
%       tan(a)            =  R(j,i) / ( -sigma * R(l,i) )
%       tan(c)            =  R(i,j) / (  sigma * R(i,l) )
% =====================================================================
function [s1, s2] = solve_proper(R, i, j, l, sigma, tol)
    cos_b = R(i,i);
    if ~isa(R, 'sym')
        cos_b = max(min(cos_b, 1), -1);         % clamp against round-off
    end
    sb_sq = R(i,j)^2 + R(i,l)^2;

    if isa(R, 'sym') || sb_sq > tol^2
        sin_b = sqrt(sb_sq);
        b1 = atan2( sin_b,  cos_b);             % b1 in (0, pi)
        a1 = atan2( R(j,i), -sigma*R(l,i));
        c1 = atan2( R(i,j),  sigma*R(i,l));

        b2 = atan2(-sin_b,  cos_b);             % b2 = -b1
        a2 = atan2(-R(j,i),  sigma*R(l,i));
        c2 = atan2(-R(i,j), -sigma*R(i,l));

        s1 = [a1, b1, c1];
        s2 = [a2, b2, c2];
    else
        % Gimbal lock: b = 0 or pi.  Same trick: residual = R_j(-b) * R
        % collapses to R_i(c), and we set a = 0.
        b1 = atan2(0, cos_b);                   % 0 or pi
        a1 = 0;
        residual = elem_rot(j, -b1) * R;
        c1 = extract_elementary(residual, i);
        s1 = [a1, b1, c1];
        s2 = s1;
    end
end


% =====================================================================
%   Small helpers
% =====================================================================
function s = perm_sign(p)
%PERM_SIGN  +1 if p is an even (cyclic) permutation of (1,2,3), -1 else.
    p = p(:).';
    if isequal(p, [1 2 3]) || isequal(p, [2 3 1]) || isequal(p, [3 1 2])
        s = 1;
    else
        s = -1;
    end
end


function R = elem_rot(axis_idx, angle)
%ELEM_ROT  Elementary rotation matrix about a principal axis.
    c = cos(angle); s = sin(angle);
    switch axis_idx
        case 1, R = [1  0  0;  0  c -s;  0  s  c];
        case 2, R = [c  0  s;  0  1  0; -s  0  c];
        case 3, R = [c -s  0;  s  c  0;  0  0  1];
    end
end


function ang = extract_elementary(R, axis_idx)
%EXTRACT_ELEMENTARY  Recover the angle from a matrix of the form R_axis(angle).
    switch axis_idx
        case 1, ang = atan2(R(3,2), R(3,3));
        case 2, ang = atan2(R(1,3), R(1,1));
        case 3, ang = atan2(R(2,1), R(1,1));
    end
end
