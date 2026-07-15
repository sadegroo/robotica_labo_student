function R = Euler2RotMat(angles, convention)
%Euler2RotMat  Build a rotation matrix from an Euler or fixed-angle triple.
%
%   R = Euler2RotMat(angles, convention) returns the 3x3 rotation
%   matrix corresponding to the given angle triple.  angles is a
%   3-element vector [theta_a, theta_b, theta_c] in radians, ordered
%   as returned by RotMat2Euler.
%
%   convention is a 4-character string:
%       characters 1-3 : axis sequence, each in {'x','y','z'}, with
%                        adjacent axes different (12 valid sequences)
%       character 4    : 'e' for Euler (rotations about the moving frame)
%                        'f' for fixed (rotations about the fixed frame)
%   Combined, this covers all 24 conventions of Craig's Appendix B.
%
%   Composition:
%       'abce' gives R = R_a(theta_a) * R_b(theta_b) * R_c(theta_c).
%       'abcf' gives R = R_c(theta_c) * R_b(theta_b) * R_a(theta_a).
%
%   This is the inverse of RotMat2Euler:
%       Euler2RotMat(RotMat2Euler(R, conv), conv) == R.
%
%   angles may be numeric or symbolic (sym); with symbolic input the
%   returned matrix is symbolic.
%
%   If angles or the convention string is invalid, R is returned as
%   NaN(3,3).

    % --- Validate angles --------------------------------------------------
    if numel(angles) ~= 3 || ~(isnumeric(angles) || isa(angles, 'sym')) ...
            || (isnumeric(angles) && ~all(isfinite(angles)))
        R = nan(3,3); return
    end

    % --- Validate convention ----------------------------------------------
    convention = lower(char(convention));
    if numel(convention) ~= 4 || ...
            ~all(ismember(convention(1:3), 'xyz')) || ...
            ~ismember(convention(4), {'e','f'}) || ...
            convention(1) == convention(2) || ...
            convention(2) == convention(3)
        R = nan(3,3); return
    end

    axes_str = convention(1:3);
    mode     = convention(4);

    % --- Map fixed -> equivalent Euler ------------------------------------
    %   'abcf' with triple (a,b,c) composes R_c * R_b * R_a, which is the
    %   Euler sequence 'cba'e with the reversed triple (c,b,a).
    if mode == 'f'
        axes_str = fliplr(axes_str);
        angles   = angles(end:-1:1);
    end

    R = elem_rot(axes_str(1), angles(1)) ...
      * elem_rot(axes_str(2), angles(2)) ...
      * elem_rot(axes_str(3), angles(3));
end


function R = elem_rot(axis_char, angle)
%ELEM_ROT  Elementary rotation matrix about a principal axis.
    c = cos(angle); s = sin(angle);
    switch axis_char
        case 'x', R = [1  0  0;  0  c -s;  0  s  c];
        case 'y', R = [c  0  s;  0  1  0; -s  0  c];
        case 'z', R = [c -s  0;  s  c  0;  0  0  1];
    end
end
