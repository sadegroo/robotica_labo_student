function ellipsoid = manipulabilityEllipsoid(J, part)
% MANIPULABILITYELLIPSOID  Manipulability ellipsoid data from a Jacobian.
%   ellipsoid = manipulabilityEllipsoid(J) takes a 6xn Jacobian J (as
%   returned by JacobianSpace, rows 1:3 angular / rows 4:6 linear) and
%   returns the data needed to plot the linear-velocity manipulability
%   ellipsoid, without plotting it.
%
%   ellipsoid = manipulabilityEllipsoid(J, part) selects the ellipsoid:
%     'linear'  - use rows 4:6 of J (default)
%     'angular' - use rows 1:3 of J
%   A 3xn Jacobian is used as-is and part is ignored.
%   A 2xn planar Jacobian is also used as-is; the result is then a 2D
%   ellipse (2x2 axes, 2x1 radii) for use with plotEllipsoid.
%
%   The ellipsoid is the image of the unit sphere in joint-velocity space,
%   {Jp*qd : ||qd|| = 1}, characterized by A = Jp*Jp'.
%
%   Output struct fields:
%     axes  - matrix whose columns are unit vectors along the
%             principal axes, longest first (right-handed set)
%     radii - semi-axis lengths, sqrt of the eigenvalues of A
%     A     - symmetric matrix Jp*Jp'
%     w     - manipulability measure prod(radii) = sqrt(det(A))

    if nargin < 2
        part = 'linear';
    end
    if size(J, 1) == 6
        switch lower(part)
            case 'linear'
                Jp = J(4:6, :);
            case 'angular'
                Jp = J(1:3, :);
            otherwise
                error('part must be ''linear'' or ''angular''.');
        end
    elseif size(J, 1) == 3 || size(J, 1) == 2
        Jp = J;
    else
        error('J must have 2, 3 or 6 rows.');
    end

    A = Jp * Jp.';
    A = (A + A.') / 2;              % symmetrize against round-off
    [V, D] = eig(A);
    lambda = max(diag(D), 0);       % clip tiny negative eigenvalues
    [lambda, order] = sort(lambda, 'descend');
    V = V(:, order);
    if det(V) < 0                   % keep the axes right-handed
        V(:, end) = -V(:, end);
    end

    ellipsoid.axes = V;
    ellipsoid.radii = sqrt(lambda);
    ellipsoid.A = A;
    ellipsoid.w = prod(ellipsoid.radii);
end
