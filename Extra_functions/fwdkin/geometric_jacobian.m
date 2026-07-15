function [Jv_geo, Jw_geo] = geometric_jacobian(DH_table, joint_types, frame_idx)
%GEOMETRIC_JACOBIAN  Column-wise closed-form Jacobian (reference for tests).
%
%   [JV, JW] = GEOMETRIC_JACOBIAN(DH_TABLE, JOINT_TYPES) returns the linear
%   and angular Jacobians (in frame {0}) for the origin of the LAST DH
%   frame, computed via the textbook column-wise formula:
%       revolute  : Jv(:,i) = z_{i-1} x (p_target - p_{i-1}),  Jw(:,i) = z_{i-1}
%       prismatic : Jv(:,i) = z_{i-1},                          Jw(:,i) = 0
%
%   [JV, JW] = GEOMETRIC_JACOBIAN(DH_TABLE, JOINT_TYPES, K) returns the same
%   for the origin of frame {K}. Columns past K are left at zero (those
%   joints come after the target frame and so cannot affect it).
%
%   The function uses the same symbolic joint variables (sym('q',[n 1]))
%   that prop_vel_jac uses, so results can be compared symbol-for-symbol
%   after substituting a numeric q.

    n = size(DH_table, 1);
    if nargin < 3 || isempty(frame_idx), frame_idx = n; end

    jt = lower(joint_types(:).');
    if numel(jt) == 1
        jt = repmat(jt, 1, n);
    end

    % Same symbolic DH table that prop_vel_jac builds
    q = sym('q', [n 1]);
    DH_sym = sym(DH_table);
    for i = 1:n
        if jt(i) == 'r'
            DH_sym(i, 4) = DH_sym(i, 4) + q(i);
        else
            DH_sym(i, 3) = DH_sym(i, 3) + q(i);
        end
    end

    [~, ~, Tcumul] = DH_full(DH_sym);
    p_target = Tcumul{frame_idx}(1:3, 4);

    Jv_geo = sym(zeros(3, n));
    Jw_geo = sym(zeros(3, n));

    for k = 1:frame_idx
        if k == 1
            p_prev = sym([0; 0; 0]);    % origin of {0}
            z_prev = sym([0; 0; 1]);    % z-axis of {0}
        else
            p_prev = Tcumul{k-1}(1:3, 4);
            z_prev = Tcumul{k-1}(1:3, 3);
        end

        if jt(k) == 'r'
            Jv_geo(:, k) = cross(z_prev, p_target - p_prev);
            Jw_geo(:, k) = z_prev;
        else
            Jv_geo(:, k) = z_prev;
            % Jw_geo(:, k) stays zero
        end
    end

    Jv_geo = simplify(Jv_geo);
    Jw_geo = simplify(Jw_geo);
end