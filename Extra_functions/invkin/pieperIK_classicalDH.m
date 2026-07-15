function solutions = pieperIK_classicalDH(DH_table, T_desired)
% PIEPERIK_CLASSICALDH  Closed-form inverse kinematics via Pieper's method
%   for a 6-DOF revolute manipulator whose last three joint axes intersect
%   at a common point. Uses CLASSICAL Denavit-Hartenberg parameters.
%
%   solutions = pieperIK_classicalDH(DH_table, T_desired)
%
%   The algorithm follows Craig, "Introduction to Robotics: Mechanics and
%   Control" (3rd ed.), Section 4.6, but the equations have been re-derived
%   from scratch because Craig uses MODIFIED DH and this routine accepts
%   CLASSICAL DH (the link transformation matrix is different, so f1, f2,
%   f3, k1..k4, and the wrist-centre position in {3} all change).
%
% INPUTS
%   DH_table  : 6x4 double, rows = [alpha_i, a_i, d_i, theta_i].
%               Angles in RADIANS. The theta_i column is ignored (those
%               are the unknowns); pass zeros there if you like.
%   T_desired : 4x4 desired homogeneous pose ^0T_6 of the end-effector.
%
% OUTPUT
%   solutions : 6 x N matrix. Each column is one joint-vector solution
%               [theta_1; theta_2; theta_3; theta_4; theta_5; theta_6].
%               N is at most 8 (4 arm configurations x 2 wrist branches).
%               Angles are wrapped to (-pi, pi].
%
% PIEPER ASSUMPTIONS (classical DH form, checked at runtime)
%   a(4)     = 0     z_3 and z_4 intersect
%   a(5)     = 0     z_4 and z_5 intersect
%   d(5)     = 0     intersection is a SINGLE common point (the wrist)
%   a(6)     = 0     tool offset is purely along z_5
%   alpha(6) = 0     z_6 = z_5  (so W = p_EE - d_6 * z_EE)
%
% CLASSICAL DH CONVENTION
%   ^{i-1}A_i = Rot(z, theta_i) * Trans(z, d_i) * Trans(x, a_i) * Rot(x, alpha_i)
%
%             | c_t   -s_t*c_a   s_t*s_a   a*c_t |
%             | s_t    c_t*c_a  -c_t*s_a   a*s_t |
%             |  0       s_a      c_a       d   |
%             |  0        0        0        1   |
%
%   This differs from Craig's modified DH, where the rotation/translation
%   in x is applied BEFORE rotation/translation in z.  As a consequence:
%
%     - the wrist centre in frame {3} is simply [0; 0; d4]
%       (Craig's modified DH gives [a3; -s_a3*d4; c_a3*d4]).
%     - frame {1} has its origin offset from {0} by (a1*c1, a1*s1, d1),
%       so the convenient base-frame radius squared is
%             R_eq = x^2 + y^2 + (z-d1)^2 - a1^2,   not  x^2+y^2+z^2.
%
% EXAMPLE (PUMA 560 - classical DH)
%   DH = [ -pi/2,    0,       0,       0;
%             0,    0.4318,   0,       0;
%           pi/2,   0.0203,   0.15005, 0;
%          -pi/2,   0,        0.4318,  0;
%           pi/2,   0,        0,       0;
%             0,    0,        0,       0 ];
%   theta_true = [10; 20; -30; 40; 50; 60]*pi/180;
%   T = fk_classicalDH(DH, theta_true);
%   sols = pieperIK_classicalDH(DH, T);
%   % verify: at least one column of sols equals theta_true (up to wrap)

    % ------------------------------------------------------------------
    % 0.  Validation
    % ------------------------------------------------------------------
    assert(isequal(size(DH_table),[6 4]), 'DH_table must be 6x4.');
    assert(isequal(size(T_desired),[4 4]), 'T_desired must be 4x4.');

    alpha = DH_table(:,1);
    a     = DH_table(:,2);
    d     = DH_table(:,3);

    tolP = 1e-9;
    if abs(a(4))     > tolP || abs(a(5))     > tolP || abs(d(5)) > tolP || ...
       abs(a(6))     > tolP || abs(alpha(6)) > tolP
        error(['Pieper''s assumptions violated. Required (classical DH): ',...
               'a(4)=a(5)=a(6)=0, d(5)=0, alpha(6)=0.']);
    end

    % ------------------------------------------------------------------
    % 1.  Wrist centre in base frame:   W = p_EE - d6 * z_EE
    % ------------------------------------------------------------------
    R_des = T_desired(1:3,1:3);
    p_des = T_desired(1:3,4);
    z_EE  = R_des(:,3);
    W     = p_des - d(6) * z_EE;
    x = W(1);  y = W(2);  z = W(3);

    % ------------------------------------------------------------------
    % 2.  Pre-computed sines, cosines, and short names
    % ------------------------------------------------------------------
    sa1 = sin(alpha(1)); ca1 = cos(alpha(1));
    sa2 = sin(alpha(2)); ca2 = cos(alpha(2));
    sa3 = sin(alpha(3)); ca3 = cos(alpha(3));
    sa4 = sin(alpha(4)); ca4 = cos(alpha(4));
    sa5 = sin(alpha(5)); ca5 = cos(alpha(5));

    a1 = a(1); a2 = a(2); a3 = a(3);
    d1 = d(1); d2 = d(2); d3 = d(3); d4 = d(4);

    f3 = d3 + d4*ca3;        % constant in theta_3

    % ------------------------------------------------------------------
    % 3.  Reduced equations   (see derivation notes)
    %
    %   With  W in {3}  =  [0 ; 0 ; d4]   (CLASSICAL DH),
    %   we get:
    %       f1 = a3 c3 + d4 sa3 s3
    %       f2 = a3 s3 - d4 sa3 c3
    %       f3 = d3 + d4 ca3                       (constant in c3,s3)
    %
    %       g1 =  c2*f1 - s2*ca2*f2 + s2*sa2*f3 + a2*c2 = c2*k1 + s2*k2
    %       g2 =  s2*f1 + c2*ca2*f2 - c2*sa2*f3 + a2*s2 = s2*k1 - c2*k2
    %       g3 =  sa2*f2 + ca2*f3 + d2                  (constant in c2,s2)
    %
    %   with k1 = a2 + f1,   k2 = sa2*f3 - ca2*f2.
    %
    %   Eliminating theta_1 (frame {1} origin is offset from {0} by
    %   (a1*c1, a1*s1, d1) under classical DH):
    %       R_eq := x^2 + y^2 + (z-d1)^2 - a1^2
    %             = g1^2 + 2*a1*g1 + g2^2 + g3^2
    %             = 2*a1*(k1*c2 + k2*s2) + k3
    %       z_eq := z - d1
    %             = sa1*(k1*s2 - k2*c2) + k4
    %
    %   k3 collapses to LINEAR in (c3,s3) thanks to f1^2 + f2^2 = const:
    %       k3 = K_const + K_c*c3 + K_s*s3
    %       k4 = ca1*(sa2*f2 + ca2*f3 + d2)
    % ------------------------------------------------------------------
    R_eq = x^2 + y^2 + (z - d1)^2 - a1^2;
    z_eq = z - d1;

    % ------------------------------------------------------------------
    % 4.  Coefficient vectors of (1+t^2)*k_i(theta_3) where t=tan(th3/2)
    %     stored as MATLAB polynomials in DESCENDING powers [t^2, t, 1]
    % ------------------------------------------------------------------
    K1 = [a2 - a3,                2*d4*sa3,         a2 + a3];
    K2 = [sa2*f3 - ca2*d4*sa3,   -2*ca2*a3,         sa2*f3 + ca2*d4*sa3];

    K_const = a3^2 + d4^2*sa3^2 + sa2^2*f3^2 + (ca2*f3 + d2)^2 + a2^2;
    K_c     = 2*(a2*a3            - sa2*d2*d4*sa3);
    K_s     = 2*(a2*d4*sa3        + sa2*d2*a3);
    K3      = [K_const - K_c,     2*K_s,            K_const + K_c];

    K4 = [ca1*( sa2*d4*sa3 + ca2*f3 + d2), ...
          2*ca1*sa2*a3, ...
          ca1*(-sa2*d4*sa3 + ca2*f3 + d2)];

    one_t2 = [1 0 1];                              % (1 + t^2)

    % U(t) := (1+t^2)*R_eq - K3(t)        (degree 2)
    % V(t) := (1+t^2)*z_eq - K4(t)        (degree 2)
    U = polyAdd(polyMul(R_eq, one_t2), -K3);
    V = polyAdd(polyMul(z_eq, one_t2), -K4);

    % ------------------------------------------------------------------
    % 5.  Solve for theta_3
    %
    %     General quartic from squaring & adding:
    %       sa1^2 * U^2  +  4 a1^2 * V^2  -  4 a1^2 sa1^2 * (K1^2 + K2^2) = 0
    % ------------------------------------------------------------------
    a1_zero  = abs(a1)  < tolP;
    sa1_zero = abs(sa1) < tolP;

    if a1_zero && sa1_zero
        warning(['Both a1 and sin(alpha1) are zero - base link is degenerate. ',...
                 'theta_1 will be set to zero; results may be incomplete.']);
        t_roots = realRoots(U);          % R_eq = k3 still required
    elseif a1_zero
        % R_eq = k3   (quadratic in t)
        t_roots = realRoots(U);
    elseif sa1_zero
        % z_eq = k4   (quadratic in t)
        t_roots = realRoots(V);
    else
        Usq         = polyMul(U, U);                 % deg 4
        Vsq         = polyMul(V, V);                 % deg 4
        K1sq_p_K2sq = polyAdd(polyMul(K1,K1), polyMul(K2,K2));   % deg 4
        P  = polyAdd(sa1^2 * Usq, 4*a1^2 * Vsq);
        P  = polyAdd(P, -4*a1^2*sa1^2 * K1sq_p_K2sq);
        t_roots = realRoots(P);
    end

    if isempty(t_roots)
        warning('Pieper IK: no real solutions for theta_3.');
        solutions = zeros(6,0);
        return;
    end

    % ------------------------------------------------------------------
    % 6.  For each theta_3 root recover theta_2 then theta_1
    % ------------------------------------------------------------------
    arm_sols = zeros(3,0);     % rows: theta_1, theta_2, theta_3

    for t = t_roots(:).'
        c3  = (1 - t^2)/(1 + t^2);
        s3  =  2*t    /(1 + t^2);
        th3 = atan2(s3, c3);

        f1 = a3*c3 + d4*sa3*s3;
        f2 = a3*s3 - d4*sa3*c3;

        k1 = a2 + f1;
        k2 = sa2*f3 - ca2*f2;
        k3 = K_const + K_c*c3 + K_s*s3;
        k4 = ca1*(sa2*f2 + ca2*f3 + d2);

        if ~a1_zero && ~sa1_zero
            % Generic case - unique (c2, s2).
            rhs1 = (R_eq - k3)/(2*a1);
            rhs2 = (z_eq - k4)/sa1;
            den  = k1^2 + k2^2;
            if den < 1e-14
                continue;
            end
            c2 = (k1*rhs1 - k2*rhs2)/den;
            s2 = (k1*rhs2 + k2*rhs1)/den;
            th2 = atan2(s2, c2);
            th1 = solveTheta1(c2, s2, c3, s3, x, y, z, alpha, a, d);
            arm_sols(:,end+1) = [th1; th2; th3]; %#ok<AGROW>

        elseif a1_zero && ~sa1_zero
            % Two (c2, s2) branches per theta_3 from
            %   -k2*c2 + k1*s2 = (z_eq - k4)/sa1   on the unit circle.
            rhs2 = (z_eq - k4)/sa1;
            [c2_a, s2_a, c2_b, s2_b, ok] = solveLinearPlusUnit(-k2, k1, rhs2);
            if ~ok, continue; end
            for branch = 1:2
                if branch == 1, c2 = c2_a; s2 = s2_a;
                else,            c2 = c2_b; s2 = s2_b; end
                th2 = atan2(s2, c2);
                th1 = solveTheta1(c2, s2, c3, s3, x, y, z, alpha, a, d);
                arm_sols(:,end+1) = [th1; th2; th3]; %#ok<AGROW>
            end

        elseif ~a1_zero && sa1_zero
            % Two (c2, s2) branches per theta_3 from
            %   k1*c2 + k2*s2 = (R_eq - k3)/(2 a1)   on the unit circle.
            rhs1 = (R_eq - k3)/(2*a1);
            [c2_a, s2_a, c2_b, s2_b, ok] = solveLinearPlusUnit(k1, k2, rhs1);
            if ~ok, continue; end
            for branch = 1:2
                if branch == 1, c2 = c2_a; s2 = s2_a;
                else,            c2 = c2_b; s2 = s2_b; end
                th2 = atan2(s2, c2);
                th1 = solveTheta1(c2, s2, c3, s3, x, y, z, alpha, a, d);
                arm_sols(:,end+1) = [th1; th2; th3]; %#ok<AGROW>
            end

        else
            % a1 = sa1 = 0  -> theta_1 underdetermined; pick zero.
            th2 = 0;  th1 = 0;
            arm_sols(:,end+1) = [th1; th2; th3]; %#ok<AGROW>
        end
    end

    if isempty(arm_sols)
        solutions = zeros(6,0);
        return;
    end

    % ------------------------------------------------------------------
    % 7.  Wrist decomposition  (theta_4, theta_5, theta_6)
    %
    %     ^3R_6 = ^0R_3' * R_des                  (call this M)
    %
    %     Comparing entries of M with the symbolic ^3R_6:
    %       M(3,3) = ca4*ca5 - sa4*sa5*c5
    %     -> c5 = (ca4*ca5 - M(3,3)) / (sa4*sa5)
    %
    %     With A = s5*sa5,                  B = ca4*c5*sa5 + sa4*ca5
    %       theta_4 = atan2( B*M(1,3) + A*M(2,3),  A*M(1,3) - B*M(2,3) )
    %
    %     With C = sa4*s5,                  D = sa4*c5*ca5 + ca4*sa5
    %       theta_6 = atan2( D*M(3,1) - C*M(3,2),  C*M(3,1) + D*M(3,2) )
    %
    %     Two branches for s5 = +/- sqrt(1 - c5^2)  give the two wrist
    %     configurations ("flip"/"no-flip") -- equation (4.1) in Craig.
    % ------------------------------------------------------------------
    if abs(sa4) < tolP || abs(sa5) < tolP
        warning(['Wrist alpha angles make the standard decomposition ',...
                 'ill-defined (alpha_4 or alpha_5 in {0, pi}). ',...
                 'A more specialised wrist solver is required.']);
    end

    out = zeros(6,0);
    warned_singular = false;

    for k = 1:size(arm_sols,2)
        th1 = arm_sols(1,k);
        th2 = arm_sols(2,k);
        th3 = arm_sols(3,k);

        R03 = dh_rot(alpha(1),th1) * dh_rot(alpha(2),th2) * dh_rot(alpha(3),th3);
        M   = R03.' * R_des;

        if abs(sa4*sa5) < tolP
            c5 = max(-1, min(1, M(3,3)));
        else
            c5 = (ca4*ca5 - M(3,3)) / (sa4*sa5);
            c5 = max(-1, min(1, c5));
        end
        s5_options = [+sqrt(max(0,1-c5^2)), -sqrt(max(0,1-c5^2))];

        for s5 = s5_options
            th5 = atan2(s5, c5);

            if abs(s5) < 1e-7
                % Wrist singularity (theta_5 ~ 0 or pi).  Set theta_4 = 0
                % and absorb the degenerate rotation into theta_6.
                if ~warned_singular
                    warning('Wrist singularity (theta_5 ~ 0 or pi); setting theta_4 = 0.');
                    warned_singular = true;
                end
                th4 = 0;
                Rtail = (dh_rot(alpha(4),th4)*dh_rot(alpha(5),th5)).' * M;
                th6   = atan2(Rtail(2,1), Rtail(1,1));
            else
                A =  s5*sa5;
                B =  ca4*c5*sa5 + sa4*ca5;
                th4 = atan2(B*M(1,3) + A*M(2,3),  A*M(1,3) - B*M(2,3));

                C =  sa4*s5;
                D =  sa4*c5*ca5 + ca4*sa5;
                th6 = atan2(D*M(3,1) - C*M(3,2),  C*M(3,1) + D*M(3,2));
            end

            out(:,end+1) = [th1; th2; th3; th4; th5; th6]; %#ok<AGROW>
        end
    end

    % Wrap to (-pi, pi].
    solutions = atan2(sin(out), cos(out));
end


% =====================================================================
% =========================  LOCAL FUNCTIONS  =========================
% =====================================================================

function R = dh_rot(alpha_i, theta_i)
% Rotation portion of the classical-DH transform ^{i-1}A_i.
    ct = cos(theta_i); st = sin(theta_i);
    ca = cos(alpha_i); sa = sin(alpha_i);
    R = [ct, -st*ca,  st*sa;
         st,  ct*ca, -ct*sa;
          0,     sa,     ca];
end


function p = polyAdd(a, b)
% Add two polynomials given in MATLAB descending-power form.
    n = max(numel(a), numel(b));
    a = [zeros(1,n-numel(a)) a];
    b = [zeros(1,n-numel(b)) b];
    p = a + b;
end


function p = polyMul(a, b)
    p = conv(a, b);
end


function r = realRoots(P)
% Real roots of polynomial P (descending powers); leading zeros trimmed.
    while numel(P) > 1 && abs(P(1)) < 1e-12
        P = P(2:end);
    end
    if numel(P) <= 1
        r = [];
        return;
    end
    rts = roots(P);
    mask = abs(imag(rts)) < 1e-6 * max(1, max(abs(real(rts))));
    r    = real(rts(mask));
end


function [c_a, s_a, c_b, s_b, ok] = solveLinearPlusUnit(A, B, C)
% Solve   A*c + B*s = C   subject to   c^2 + s^2 = 1.
% Returns the two solutions (c_a,s_a) and (c_b,s_b).  ok = false if none.
    R2 = A^2 + B^2;
    if R2 < 1e-14
        ok = false; c_a=0; s_a=0; c_b=0; s_b=0; return;
    end
    R = sqrt(R2);
    if abs(C) > R + 1e-9
        ok = false; c_a=0; s_a=0; c_b=0; s_b=0; return;
    end
    ok = true;
    beta  = atan2(B, A);
    delta = acos(max(-1, min(1, C/R)));
    g_a = beta + delta;     c_a = cos(g_a);  s_a = sin(g_a);
    g_b = beta - delta;     c_b = cos(g_b);  s_b = sin(g_b);
end


function th1 = solveTheta1(c2, s2, c3, s3, x, y, z, alpha, a, d)
% Recover theta_1 once theta_2 and theta_3 are known.
%
%   x = c1*A1 + s1*B1,   y = s1*A1 - c1*B1
%   with  A1 = g1 + a1,   B1 = -ca1*g2 + sa1*g3
    sa1 = sin(alpha(1)); ca1 = cos(alpha(1));
    sa2 = sin(alpha(2)); ca2 = cos(alpha(2));
    sa3 = sin(alpha(3)); ca3 = cos(alpha(3));
    a1  = a(1); a2 = a(2); a3 = a(3);
    d2  = d(2); d3 = d(3); d4 = d(4);

    f1 = a3*c3 + d4*sa3*s3;
    f2 = a3*s3 - d4*sa3*c3;
    f3 = d3 + d4*ca3;

    g1 =  c2*f1 - s2*ca2*f2 + s2*sa2*f3 + a2*c2;
    g2 =  s2*f1 + c2*ca2*f2 - c2*sa2*f3 + a2*s2;
    g3 =  sa2*f2 + ca2*f3 + d2;

    A1 =  g1 + a1;
    B1 = -ca1*g2 + sa1*g3;
    den = A1^2 + B1^2;

    if den < 1e-14
        th1 = 0;     % wrist centre on the z_0 axis -- theta_1 free.
        return;
    end
    c1 = (A1*x - B1*y)/den;
    s1 = (B1*x + A1*y)/den;
    th1 = atan2(s1, c1);
end


function T = fk_classicalDH(DH_table, theta)
% Forward kinematics with classical DH (helper for verification).
% theta: nx1 vector of joint angles (radians). Overrides DH_table(:,4).
    T = eye(4);
    for i = 1:size(DH_table,1)
        al = DH_table(i,1);
        ai = DH_table(i,2);
        di = DH_table(i,3);
        th = theta(i);
        ct = cos(th); st = sin(th);
        ca = cos(al); sa = sin(al);
        Ai = [ct, -st*ca,  st*sa,  ai*ct;
              st,  ct*ca, -ct*sa,  ai*st;
               0,     sa,     ca,     di;
               0,      0,      0,      1];
        T = T * Ai;
    end
end