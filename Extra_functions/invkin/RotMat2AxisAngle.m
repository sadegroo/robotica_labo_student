function [K1, theta1, K2, theta2] = RotMat2AxisAngle(R)
%ROTMAT_TO_AXIS_ANGLE  Extract the two equivalent angle-axis representations
%   of a 3x3 rotation matrix.
%
%   [K1, theta1, K2, theta2] = rotmat_to_axis_angle(R) returns:
%       K1, theta1 : a unit axis (3x1) and angle (radians) such that
%                    R = R_K1(theta1)
%       K2, theta2 : the equivalent pair (-K1, -theta1) describing the
%                    same rotation
%
%   If R is not a proper orthonormal matrix (R'*R = I and det(R) = +1),
%   all outputs are set to NaN.
%
%   Special cases:
%     theta = 0 : R = I, no rotation. Axis is undefined; K1 returned as NaN.
%     theta = pi: sin(theta) = 0, antisymmetric formula fails. K is
%                 recovered from the symmetric part R + I = 2*K*K'.
%
%   Reference: Craig, "Introduction to Robotics", eq. (2.80) and (2.82).

    tol = 1e-6;

    % ---- 1. Validate input ------------------------------------------------
    if ~isequal(size(R), [3 3]) ...
            || norm(R.' * R - eye(3), 'fro') > tol ...
            || abs(det(R) - 1) > tol
        K1 = nan(3,1);  theta1 = NaN;
        K2 = nan(3,1);  theta2 = NaN;
        return
    end

    % ---- 2. Recover theta from the trace ---------------------------------
    %   trace(R) = 1 + 2*cos(theta)
    cos_theta = (trace(R) - 1) / 2;
    cos_theta = max(min(cos_theta, 1), -1);   % clamp against round-off
    theta1    = acos(cos_theta);              % in [0, pi]

    % ---- 3. Recover the unit axis K --------------------------------------
    if sin(theta1) > tol
        % Generic case: use the antisymmetric part
        %   R - R' = 2*sin(theta) * [K]_x
        K1 = (1 / (2 * sin(theta1))) * [ R(3,2) - R(2,3);
                                          R(1,3) - R(3,1);
                                          R(2,1) - R(1,2) ];

    elseif theta1 < tol
        % theta = 0 : identity rotation, axis is arbitrary
        K1 = nan(3,1);

    else
        % theta = pi : antisymmetric part vanishes, use symmetric part
        %   R + I = 2*K*K'   so   M(i,j) = k_i * k_j
        M = (R + eye(3)) / 2;

        % Pick the largest diagonal entry to compute one component robustly
        diag_vals = max(diag(M), 0);          % clamp against round-off
        [~, idx]  = max(diag_vals);

        K1      = zeros(3,1);
        K1(idx) = sqrt(diag_vals(idx));

        % Other two components from the off-diagonal entries of M
        for i = 1:3
            if i ~= idx
                K1(i) = M(idx, i) / K1(idx);
            end
        end

        % Re-normalise to kill any residual numerical drift
        K1 = K1 / norm(K1);
    end

    % ---- 4. Second equivalent solution -----------------------------------
    %   (K, theta) and (-K, -theta) describe the same rotation matrix.
    K2     = -K1;
    theta2 = -theta1;
end