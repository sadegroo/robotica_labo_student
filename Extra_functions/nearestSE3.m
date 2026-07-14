function [T, d] = nearestSE3(Ttilde)
% NEARESTSE3  Project Ttilde onto SE(3) in Frobenius norm.
%   [T,d] = nearestSE3(Ttilde) returns the closest homogeneous transform T
%   in SE(3) to the 4x4 matrix Ttilde and the Frobenius distance
%   d = ||Ttilde - T||_F.
%
%   By the block decomposition (see derivation, Step 13):
%       d^2 = ||M - R*||_F^2  +  ||p - t*||^2  +  ||lastRow - [0 0 0 1]||^2
%           = ||M - R*||_F^2  +  0             +  ||lastRow - [0 0 0 1]||^2
%   since t* = p is unconstrained.

    M  = Ttilde(1:3, 1:3);
    p  = Ttilde(1:3, 4);
    lr = Ttilde(4, :);                      % last row of Ttilde

    [R, dR] = nearestSO3(M);

    T = [R, p; 0 0 0 1];

    dLast = norm(lr - [0 0 0 1]);           % penalty if bottom row drifted
    d = sqrt(dR^2 + dLast^2);
end