function [R, d] = nearestSO3(M)
% NEARESTSO3  Project M onto SO(3) in Frobenius norm.
%   [R,d] = nearestSO3(M) returns the closest rotation R in SO(3) to the
%   3x3 matrix M and the Frobenius distance d = ||M - R||_F.

    [U, ~, V] = svd(M);
    D = eye(3);
    D(3,3) = det(U * V.');   % +/- 1
    R = U * D * V.';
    d = norm(M - R, 'fro');
end