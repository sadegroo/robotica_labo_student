function R = RotMatKaxis(k, theta)
% RotMatKaxis  Rotation matrix from angle-axis (Rodrigues') formula.
%   R = RotMatKaxis(k, theta) returns the 3x3 rotation matrix that
%   rotates by angle theta (radians) about the axis given by 3-vector k.
%   The axis k does not need to be unit length; it is normalised here.
%
%   Implements Craig (2.80):
%       R = I*c + (1-c)*k*k' + s*[k]_x
%   with c = cos(theta), s = sin(theta), v = 1 - c (versine),
%   and [k]_x the skew-symmetric matrix of the unit axis.

    % --- input checks & normalise axis ----------------------------------
    k = k(:);
    if numel(k) ~= 3
        error('RotMatZaxis:badAxis', 'Axis k must be a 3-vector.');
    end
    nk = norm(k);
    if nk < eps
        error('RotMatZaxis:zeroAxis', 'Axis k must be non-zero.');
    end
    k = k / nk;

    % --- shorthand ------------------------------------------------------
    c  = cos(theta);
    s  = sin(theta);
    v  = 1 - c;                       % versine
    kx = k(1);  ky = k(2);  kz = k(3);

    % --- Craig eq. (2.80), written out element-by-element ---------------
    R = [ kx*kx*v + c    ,  kx*ky*v - kz*s ,  kx*kz*v + ky*s ;
          kx*ky*v + kz*s ,  ky*ky*v + c    ,  ky*kz*v - kx*s ;
          kx*kz*v - ky*s ,  ky*kz*v + kx*s ,  kz*kz*v + c    ];
end