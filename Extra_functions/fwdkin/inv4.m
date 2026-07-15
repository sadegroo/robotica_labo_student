function [Tinv] = inv4(T)
%INV4 Compute the inverse of a 4x4 homogeneous transformation matrix.
%   Tinv = inv4(T) returns the inverse using R' and -R'*d,
%   which is more efficient and numerically stable than inv(T).

assert(isequal(size(T), [4 4]), 'Input must be a 4x4 matrix.');
assert(norm(T(4,:) - [0 0 0 1]) < 1e-10, 'Last row must be [0 0 0 1].');

R = T(1:3, 1:3);
assert(norm(R' * R - eye(3)) < 1e-10, 'Rotation part must be orthonormal.');

d = T(1:3, 4);

Tinv = [R'  -R'*d;
        0 0 0  1];

end
