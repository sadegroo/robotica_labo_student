function solutions = pieperIK_modifiedDH(DH_table, T_desired)
% PIEPERIK_MODIFIEDDH  Closed-form inverse kinematics via Pieper's method
%   for a 6-DOF revolute manipulator whose last three joint axes intersect
%   at a common point. Uses MODIFIED Denavit-Hartenberg parameters.
%
%   solutions = pieperIK_modifiedDH(DH_table, T_desired)
%
%   The algorithm follows Craig, "Introduction to Robotics: Mechanics and
%   Control" (3rd ed.), Section 4.6, which is stated directly in MODIFIED
%   DH.  In particular the reduced position equations use Craig's
%   f1..f3 (eq. 4.65) and k1..k4 (eq. 4.70).  See pieperIK_classicalDH
%   for the classical-DH counterpart.
%
% INPUTS
%   DH_table  : 6x4 double, rows = [alpha_(i-1), a_(i-1), d_i, theta_i]
%               (the same column convention as DH_modified).
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
% PIEPER ASSUMPTIONS (modified DH form, checked at runtime)
%   a_4 = 0   i.e. DH_table(5,2)     z_4 and z_5 intersect
%   a_5 = 0   i.e. DH_table(6,2)     z_5 and z_6 intersect
%   d_5 = 0   i.e. DH_table(5,3)     intersection is a SINGLE common
%                                    point (the wrist)
%   A tool offset d_6 along z_6 is allowed (W = p_EE - d_6 * z_EE);
%   modified DH has no alpha_6 or a_6 parameters in a 6-row table.
%   A non-trivial base screw (alpha_0, a_0) is allowed: the target pose
%   is pre-multiplied by [Rot(x,alpha_0)*Trans(x,a_0)]^-1 first.
%
% MODIFIED DH CONVENTION
%   ^{i-1}T_i = Rot(x, alpha_(i-1)) * Trans(x, a_(i-1))
%             * Rot(z, theta_i)     * Trans(z, d_i)
%
%   The wrist centre is the common origin of frames {4}, {5}, {6}:
%       ^3P_4ORG = [a3; -sa3*d4; ca3*d4]      (Craig eq. 4.64)
%   which yields (Craig eq. 4.65, with f3bar = d3 + d4*ca3 folded in):
%       f1 = a3 c3 + d4 sa3 s3 + a2
%       f2 = ca2*(a3 s3 - d4 sa3 c3) - sa2*f3bar
%       f3 = sa2*(a3 s3 - d4 sa3 c3) + ca2*f3bar
%   and (Craig eq. 4.70):
%       k1 = f1
%       k2 = -f2
%       k3 = f1^2 + f2^2 + f3^2 + a1^2 + d2^2 + 2 d2 f3   (linear in c3,s3)
%       k4 = ca1*(f3 + d2)
%   so that, with r^2 = x^2 + y^2 + (z-d1)^2 of the wrist centre,
%       r^2     = 2 a1 (k1 c2 + k2 s2) + k3
%       z - d1  = sa1 (k1 s2 - k2 c2) + k4
%
% EXAMPLE (PUMA 560 - modified DH, Craig ch. 3)
%   DH = [    0,   0,       0,       0;
%         -pi/2,   0,       0,       0;
%             0,   0.4318,  0.15005, 0;
%         -pi/2,   0.0203,  0.4318,  0;
%          pi/2,   0,       0,       0;
%         -pi/2,   0,       0,       0 ];
%   theta_true = [10; 20; -30; 40; 50; 60]*pi/180;
%   tbl = DH; tbl(:,4) = theta_true;
%   T = DH_modified(tbl);
%   sols = pieperIK_modifiedDH(DH, T);
%   % verify: at least one column of sols equals theta_true (up to wrap)

    % ------------------------------------------------------------------
    % 0.  Validation
    % ------------------------------------------------------------------
    assert(isequal(size(DH_table),[6 4]), 'DH_table must be 6x4.');
    assert(isequal(size(T_desired),[4 4]), 'T_desired must be 4x4.');

    % Row i holds [alpha_(i-1), a_(i-1), d_i]; index by SUBSCRIPT here:
    %   alpha(i) = alpha_(i-1),  a(i) = a_(i-1),  d(i) = d_i.
    alpha_prev = DH_table(:,1);
    a_prev     = DH_table(:,2);
    d          = DH_table(:,3);

    tolP = 1e-9;
    if abs(a_prev(5)) > tolP || abs(a_prev(6)) > tolP || abs(d(5)) > tolP
        error(['Pieper''s assumptions violated. Required (modified DH): ',...
               'a_4 = DH_table(5,2) = 0, a_5 = DH_table(6,2) = 0, ',...
               'd_5 = DH_table(5,3) = 0.']);
    end

    % ------------------------------------------------------------------
    % 1.  Remove the base screw and find the wrist centre
    %       T_eff = [Rot(x,alpha_0)*Trans(x,a_0)]^-1 * T_desired
    %       W     = p_eff - d6 * z_eff        (common origin of {4},{5},{6})
    % ------------------------------------------------------------------
    Bx = xScrew(alpha_prev(1), a_prev(1));
    T_eff = Bx \ T_desired;

    R_des = T_eff(1:3,1:3);
    p_des = T_eff(1:3,4);
    z_EE  = R_des(:,3);
    W     = p_des - d(6) * z_EE;
    x = W(1);  y = W(2);  z = W(3);

    % ------------------------------------------------------------------
    % 2.  Pre-computed sines, cosines, and short names
    %     (a1 here means Craig's a_1 = DH_table(2,2), etc.)
    % ------------------------------------------------------------------
    sa1 = sin(alpha_prev(2)); ca1 = cos(alpha_prev(2));
    sa2 = sin(alpha_prev(3)); ca2 = cos(alpha_prev(3));
    sa3 = sin(alpha_prev(4)); ca3 = cos(alpha_prev(4));
    sa4 = sin(alpha_prev(5)); ca4 = cos(alpha_prev(5));
    sa5 = sin(alpha_prev(6)); ca5 = cos(alpha_prev(6));

    a1 = a_prev(2); a2 = a_prev(3); a3 = a_prev(4);
    d1 = d(1); d2 = d(2); d3 = d(3); d4 = d(4);

    f3bar = d3 + d4*ca3;     % constant part of f3 (and of -f2/sa2)

    % ------------------------------------------------------------------
    % 3.  Reduced equations (Craig 4.64-4.70, see header)
    %       r_eq := x^2 + y^2 + (z-d1)^2 = 2 a1 (k1 c2 + k2 s2) + k3
    %       z_eq := z - d1               = sa1 (k1 s2 - k2 c2) + k4
    %     k3 is LINEAR in (c3, s3):  k3 = K_const + K_c c3 + K_s s3.
    % ------------------------------------------------------------------
    r_eq = x^2 + y^2 + (z - d1)^2;
    z_eq = z - d1;

    % ------------------------------------------------------------------
    % 4.  Coefficient vectors of (1+t^2)*k_i(theta_3) where t=tan(th3/2)
    %     stored as MATLAB polynomials in DESCENDING powers [t^2, t, 1]
    % ------------------------------------------------------------------
    K1 = [a2 - a3,               2*d4*sa3,      a2 + a3];
    K2 = [sa2*f3bar - ca2*d4*sa3, -2*ca2*a3,    sa2*f3bar + ca2*d4*sa3];

    K_const = a1^2 + a2^2 + a3^2 + d4^2*sa3^2 + f3bar^2 + d2^2 ...
              + 2*d2*ca2*f3bar;
    K_c     = 2*(a2*a3       - sa2*d2*d4*sa3);
    K_s     = 2*(a2*d4*sa3   + sa2*d2*a3);
    K3      = [K_const - K_c,    2*K_s,         K_const + K_c];

    K4 = [ca1*( sa2*d4*sa3 + ca2*f3bar + d2), ...
          2*ca1*sa2*a3, ...
          ca1*(-sa2*d4*sa3 + ca2*f3bar + d2)];

    one_t2 = [1 0 1];                              % (1 + t^2)

    % U(t) := (1+t^2)*r_eq - K3(t)        (degree 2)
    % V(t) := (1+t^2)*z_eq - K4(t)        (degree 2)
    U = polyAdd(polyMul(r_eq, one_t2), -K3);
    V = polyAdd(polyMul(z_eq, one_t2), -K4);

    % ------------------------------------------------------------------
    % 5.  Solve for theta_3
    %
    %     General quartic from squaring & adding (Craig 4.73):
    %       sa1^2 * U^2  +  4 a1^2 * V^2  -  4 a1^2 sa1^2 * (K1^2 + K2^2) = 0
    % ------------------------------------------------------------------
    a1_zero  = abs(a1)  < tolP;
    sa1_zero = abs(sa1) < tolP;

    if a1_zero && sa1_zero
        warning(['Both a1 and sin(alpha1) are zero - base link is degenerate. ',...
                 'theta_1 will be set to zero; results may be incomplete.']);
        t_roots = realRoots(U);          % r_eq = k3 still required
    elseif a1_zero
        % r_eq = k3   (quadratic in t, Craig 4.71)
        t_roots = realRoots(U);
    elseif sa1_zero
        % z_eq = k4   (quadratic in t, Craig 4.72)
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

        [f1, f2, f3] = craigF(c3, s3, a2, a3, d4, sa2, ca2, sa3, f3bar);

        k1 =  f1;
        k2 = -f2;
        k3 = K_const + K_c*c3 + K_s*s3;
        k4 = ca1*(f3 + d2);

        if ~a1_zero && ~sa1_zero
            % Generic case - unique (c2, s2).
            rhs1 = (r_eq - k3)/(2*a1);
            rhs2 = (z_eq - k4)/sa1;
            den  = k1^2 + k2^2;
            if den < 1e-14
                continue;
            end
            c2 = (k1*rhs1 - k2*rhs2)/den;
            s2 = (k2*rhs1 + k1*rhs2)/den;
            th2 = atan2(s2, c2);
            th1 = solveTheta1(c2, s2, f1, f2, f3, x, y, a1, d2, sa1, ca1);
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
                th1 = solveTheta1(c2, s2, f1, f2, f3, x, y, a1, d2, sa1, ca1);
                arm_sols(:,end+1) = [th1; th2; th3]; %#ok<AGROW>
            end

        elseif ~a1_zero && sa1_zero
            % Two (c2, s2) branches per theta_3 from
            %   k1*c2 + k2*s2 = (r_eq - k3)/(2 a1)   on the unit circle.
            rhs1 = (r_eq - k3)/(2*a1);
            [c2_a, s2_a, c2_b, s2_b, ok] = solveLinearPlusUnit(k1, k2, rhs1);
            if ~ok, continue; end
            for branch = 1:2
                if branch == 1, c2 = c2_a; s2 = s2_a;
                else,            c2 = c2_b; s2 = s2_b; end
                th2 = atan2(s2, c2);
                th1 = solveTheta1(c2, s2, f1, f2, f3, x, y, a1, d2, sa1, ca1);
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
    %     ^3R_6 = ^0R_3' * R_des,  and peeling off the leading x-screw:
    %       N := Rot(x,alpha_3)' * ^3R_6
    %          = Rz(th4) Rx(alpha_4) Rz(th5) Rx(alpha_5) Rz(th6)
    %
    %     which is the same ZXZXZ structure as in the classical solver:
    %       N(3,3) = ca4*ca5 - sa4*sa5*c5
    %     -> c5 = (ca4*ca5 - N(3,3)) / (sa4*sa5)
    %
    %     With A = s5*sa5,                  B = ca4*c5*sa5 + sa4*ca5
    %       theta_4 = atan2( B*N(1,3) + A*N(2,3),  A*N(1,3) - B*N(2,3) )
    %
    %     With C = sa4*s5,                  D = sa4*c5*ca5 + ca4*sa5
    %       theta_6 = atan2( D*N(3,1) - C*N(3,2),  C*N(3,1) + D*N(3,2) )
    %
    %     Two branches for s5 = +/- sqrt(1 - c5^2)  give the two wrist
    %     configurations ("flip"/"no-flip") -- equation (4.1) in Craig.
    % ------------------------------------------------------------------
    if abs(sa4) < tolP || abs(sa5) < tolP
        warning(['Wrist alpha angles make the standard decomposition ',...
                 'ill-defined (alpha_4 or alpha_5 in {0, pi}). ',...
                 'A more specialised wrist solver is required.']);
    end

    Rx3 = xRot(alpha_prev(4));

    out = zeros(6,0);
    warned_singular = false;

    for k = 1:size(arm_sols,2)
        th1 = arm_sols(1,k);
        th2 = arm_sols(2,k);
        th3 = arm_sols(3,k);

        % ^0R_3 = Rz(th1) * [Rx(a1) Rz(th2)] * [Rx(a2) Rz(th3)]
        R03 = dh_rot(0,th1) * dh_rot(alpha_prev(2),th2) * dh_rot(alpha_prev(3),th3);
        N   = Rx3.' * (R03.' * R_des);

        if abs(sa4*sa5) < tolP
            c5 = max(-1, min(1, N(3,3)));
        else
            c5 = (ca4*ca5 - N(3,3)) / (sa4*sa5);
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
                Rtail = (dh_rot(0,th4) * xRot(alpha_prev(5)) * ...
                         dh_rot(0,th5) * xRot(alpha_prev(6))).' * N;
                th6   = atan2(Rtail(2,1), Rtail(1,1));
            else
                A =  s5*sa5;
                B =  ca4*c5*sa5 + sa4*ca5;
                th4 = atan2(B*N(1,3) + A*N(2,3),  A*N(1,3) - B*N(2,3));

                C =  sa4*s5;
                D =  sa4*c5*ca5 + ca4*sa5;
                th6 = atan2(D*N(3,1) - C*N(3,2),  C*N(3,1) + D*N(3,2));
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

function [f1, f2, f3] = craigF(c3, s3, a2, a3, d4, sa2, ca2, sa3, f3bar)
% Craig's f1..f3 (eq. 4.65): the wrist centre seen from frame {1},
% expressed as functions of theta_3 only.
    u  = a3*c3 + d4*sa3*s3;          % x-component before the a2 shift
    v  = a3*s3 - d4*sa3*c3;
    f1 = u + a2;
    f2 = ca2*v - sa2*f3bar;
    f3 = sa2*v + ca2*f3bar;
end


function R = dh_rot(alpha_prev_i, theta_i)
% Rotation portion of the modified-DH transform ^{i-1}T_i:
%   ^{i-1}R_i = Rot(x, alpha_(i-1)) * Rot(z, theta_i)
    ct = cos(theta_i); st = sin(theta_i);
    ca = cos(alpha_prev_i); sa = sin(alpha_prev_i);
    R = [   ct,   -st,    0;
         st*ca, ct*ca,  -sa;
         st*sa, ct*sa,   ca];
end


function R = xRot(alpha)
% Elementary rotation about x.
    ca = cos(alpha); sa = sin(alpha);
    R = [1  0   0;
         0 ca -sa;
         0 sa  ca];
end


function T = xScrew(alpha, a)
% Rot(x, alpha) * Trans(x, a) as a homogeneous transform.
    T = [xRot(alpha), [a; 0; 0]; 0 0 0 1];
end


function th1 = solveTheta1(c2, s2, f1, f2, f3, x, y, a1, d2, sa1, ca1)
% Recover theta_1 once theta_2 and theta_3 are known (Craig 4.66-4.68):
%   x = c1*g1 - s1*g2,   y = s1*g1 + c1*g2
%   g1 = c2*f1 - s2*f2 + a1
%   g2 = ca1*(s2*f1 + c2*f2) - sa1*(f3 + d2)
    g1 = c2*f1 - s2*f2 + a1;
    g2 = ca1*(s2*f1 + c2*f2) - sa1*(f3 + d2);
    den = g1^2 + g2^2;

    if den < 1e-14
        th1 = 0;     % wrist centre on the z_1 axis -- theta_1 free.
        return;
    end
    c1 = (g1*x + g2*y)/den;
    s1 = (g1*y - g2*x)/den;
    th1 = atan2(s1, c1);
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
