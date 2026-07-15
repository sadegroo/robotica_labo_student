function [linvel, angvel, Jv, Jw, symbols] = prop_vel_jac(DH_table, joint_types, varargin)
%PROP_VEL_JAC  Classical-DH velocity propagation and per-link Jacobians.
%
%   [LINVEL, ANGVEL, JV, JW] = PROP_VEL_JAC(DH_TABLE, JOINT_TYPES) propagates
%   linear and angular velocities through a serial manipulator described
%   by DH_TABLE (n-by-4 double, columns [alpha a d theta] in radians) and
%   JOINT_TYPES (char vector of 'R'/'P', case-insensitive; a single 'R' or
%   'P' is broadcast to every joint).
%
%   The d and theta entries of DH_TABLE are interpreted as FIXED OFFSETS
%   on the joint variable q_i:
%       revolute  : theta_i = DH_TABLE(i,4) + q_i
%       prismatic : d_i     = DH_TABLE(i,3) + q_i
%
%   Optional name-value arguments (parsed with inputParser):
%       'q'    n-vector of joint positions  -- triggers numeric subs
%       'qdot' n-vector of joint rates      -- triggers numeric subs
%       'ee'   3-vector locating an end-effector point in the last frame
%   The literal string 'sym' may be passed as the first optional arg to
%   make symbolic intent explicit; it is the default and otherwise ignored.
%
%   Outputs:
%       LINVEL  struct with .own and .zero cell arrays; cell i is the
%               velocity of the origin of frame {i} expressed in {i} and
%               in {0} respectively. If 'ee' was given, an extra entry is
%               appended at the end of both cells.
%       ANGVEL  same layout for angular velocities (EE entry equals the
%               last link's angular velocity).
%       JV      struct with .own and .zero cell arrays; cell i is the
%               linear-velocity Jacobian for the origin of frame {i},
%               expressed in {i} or {0} respectively. EE Jacobian is
%               appended if 'ee' was provided.
%       JW      same layout for the angular-velocity Jacobians.
%   SYMBOLS     struct with q and qdot fields containing the symbolic 
%               arguments for later use in e.g. subs()
%
%   This function relies on DH_full to build the homogeneous transforms.

    % ------------------------------------------------------------------
    % Required-argument validation
    % ------------------------------------------------------------------
    if ~isnumeric(DH_table) || size(DH_table, 2) ~= 4
        error('DH_table must be an n-by-4 numeric matrix.');
    end
    n = size(DH_table, 1);

    if ~ischar(joint_types)
        error('joint_types must be a char vector.');
    end
    jt = lower(joint_types(:).');
    if numel(jt) == 1
        jt = repmat(jt, 1, n);
    elseif numel(jt) ~= n
        error('joint_types must have length 1 or %d (got %d).', n, numel(jt));
    end
    if ~all(jt == 'r' | jt == 'p')
        error('joint_types must contain only R/r and P/p characters.');
    end

    % ------------------------------------------------------------------
    % Optional-argument parsing
    % ------------------------------------------------------------------
    % Allow a bare 'sym' as the first optional arg before name-value pairs.
    if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1}, 'sym')
        varargin(1) = [];
    end

    p = inputParser;
    p.addParameter('q',    [], @(x) isempty(x) || isnumeric(x));
    p.addParameter('qdot', [], @(x) isempty(x) || isnumeric(x));
    p.addParameter('ee',   [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
    p.parse(varargin{:});

    qvals    = p.Results.q;
    qdotvals = p.Results.qdot;
    ee_pt    = p.Results.ee;
    do_subs  = ~isempty(qvals) && ~isempty(qdotvals);

    if do_subs && (numel(qvals) ~= n || numel(qdotvals) ~= n)
        error('q and qdot must each have length %d.', n);
    end

    % ------------------------------------------------------------------
    % Build symbolic joint variables and a symbolic DH table
    % ------------------------------------------------------------------
    q    = sym('q',    [n 1]);
    qdot = sym('qdot', [n 1]);

    symbols.q = q;
    symbols.qdot = qdot;

    DH_sym = sym(DH_table);
    for i = 1:n
        if jt(i) == 'r'
            DH_sym(i, 4) = DH_sym(i, 4) + q(i);
        else
            DH_sym(i, 3) = DH_sym(i, 3) + q(i);
        end
    end

    % ------------------------------------------------------------------
    % Pull rotation matrices and origin offsets from DH_full
    % ------------------------------------------------------------------
    [~, Tparts, Tcumul] = DH_full(DH_sym);
    R = cell(1, n);    % R{i} = ^{i-1}_i R
    P = cell(1, n);    % P{i} = ^{i-1} P_{i,org}
    for i = 1:n
        R{i} = Tparts{i}(1:3, 1:3);
        P{i} = Tparts{i}(1:3, 4);
    end

    % ------------------------------------------------------------------
    % Velocity propagation through the chain
    % ------------------------------------------------------------------
    w_own = cell(1, n);
    v_own = cell(1, n);

    z      = [0; 0; 1];
    w_prev = sym(zeros(3, 1));   % ^0 omega_0 = 0
    v_prev = sym(zeros(3, 1));   % ^0 v_0     = 0

    for i = 1:n
        % omega of frame {i}, still expressed in frame {i-1}
        if jt(i) == 'r'
            w_im = w_prev + qdot(i) * z;
        else
            w_im = w_prev;
        end

        % rigid-body part of the linear velocity, still in frame {i-1}
        v_im = v_prev + cross(w_im, P{i});
        if jt(i) == 'p'
            v_im = v_im + qdot(i) * z;
        end

        % rotate both into frame {i}
        w_own{i} = simplify(R{i}.' * w_im);
        v_own{i} = simplify(R{i}.' * v_im);

        w_prev = w_own{i};
        v_prev = v_own{i};
    end

    % ------------------------------------------------------------------
    % Express the same velocities in frame {0}
    % ------------------------------------------------------------------
    w_zero = cell(1, n);
    v_zero = cell(1, n);
    for i = 1:n
        R0i       = Tcumul{i}(1:3, 1:3);
        w_zero{i} = simplify(R0i * w_own{i});
        v_zero{i} = simplify(R0i * v_own{i});
    end

    % ------------------------------------------------------------------
    % Optional end-effector point (fixed in the last frame)
    % ------------------------------------------------------------------
    has_ee = ~isempty(ee_pt);
    if has_ee
        Q   = sym(ee_pt(:));
        R0n = Tcumul{n}(1:3, 1:3);

        v_ee_own  = simplify(v_own{n} + cross(w_own{n}, Q));   % image 4
        w_ee_own  = w_own{n};                                  % unchanged
        v_ee_zero = simplify(R0n * v_ee_own);
        w_ee_zero = w_zero{n};

        v_own{end+1}  = v_ee_own;
        v_zero{end+1} = v_ee_zero;
        w_own{end+1}  = w_ee_own;
        w_zero{end+1} = w_ee_zero;
    end

    % ------------------------------------------------------------------
    % Per-frame Jacobians (linear & angular, in own frame and in {0})
    % Computed for every entry in the velocity arrays, so the EE entry
    % is included automatically if it was appended above.
    % ------------------------------------------------------------------
    N       = numel(v_own);
    Jv_own  = cell(1, N);
    Jv_zero = cell(1, N);
    Jw_own  = cell(1, N);
    Jw_zero = cell(1, N);

    for i = 1:N
        Jv_own{i}  = simplify(jacobian(v_own{i},  qdot));
        Jv_zero{i} = simplify(jacobian(v_zero{i}, qdot));
        Jw_own{i}  = simplify(jacobian(w_own{i},  qdot));
        Jw_zero{i} = simplify(jacobian(w_zero{i}, qdot));
    end

    % ------------------------------------------------------------------
    % Substitute numerical q / qdot if requested
    % ------------------------------------------------------------------
    if do_subs
        old = [q; qdot];
        new = [qvals(:); qdotvals(:)];
        for i = 1:N
            v_own{i}   = double(subs(v_own{i},   old, new));
            v_zero{i}  = double(subs(v_zero{i},  old, new));
            w_own{i}   = double(subs(w_own{i},   old, new));
            w_zero{i}  = double(subs(w_zero{i},  old, new));
            Jv_own{i}  = double(subs(Jv_own{i},  old, new));
            Jv_zero{i} = double(subs(Jv_zero{i}, old, new));
            Jw_own{i}  = double(subs(Jw_own{i},  old, new));
            Jw_zero{i} = double(subs(Jw_zero{i}, old, new));
        end
    end

    % ------------------------------------------------------------------
    % Assemble output structs
    % ------------------------------------------------------------------
    linvel = struct('own', {v_own},  'zero', {v_zero});
    angvel = struct('own', {w_own},  'zero', {w_zero});
    Jv     = struct('own', {Jv_own}, 'zero', {Jv_zero});
    Jw     = struct('own', {Jw_own}, 'zero', {Jw_zero});
end