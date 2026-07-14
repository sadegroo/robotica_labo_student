function [xvec,yvec,zvec,origin] = plotFrame(MATRIX,origin)
%Plots a frame of a given input matrix and origin.
%MATRIX is either a 3x3 rotation matrix (origin given separately, defaults
%to [0;0;0]) or a 4x4 homogeneous transformation matrix (origin input is
%ignored; rotation and origin are taken from the matrix itself).
%Gives x-y-z-vector and origin as output.

tol = 1e-6;

if isequal(size(MATRIX),[4 4])
    if nargin > 1
        warning('plotFrame:originIgnored', ...
            'MATRIX is 4x4: origin input is ignored, translation part of MATRIX is used.');
    end
    R = MATRIX(1:3,1:3);
    origin = MATRIX(1:3,4);
elseif isequal(size(MATRIX),[3 3])
    R = MATRIX;
    if nargin < 2
        origin = [0;0;0];
    end
else
    error('plotFrame:badSize','MATRIX must be 3x3 or 4x4.');
end

if norm(R'*R - eye(3)) > tol
    warning('plotFrame:notOrthonormal', ...
        'Rotation part of MATRIX is not orthonormal within tolerance %g.',tol);
elseif det(R) < 0
    warning('plotFrame:leftHanded', ...
        'Rotation part of MATRIX is orthonormal but left-handed (det = %g).',det(R));
end

xvec = R(:,1);
yvec = R(:,2);
zvec = R(:,3);

%figure;
grid on
hold on
arrow3D(origin,xvec,'r');
hold on
arrow3D(origin,yvec,'g');
hold on
arrow3D(origin,zvec,'b');
xlabel('x'),ylabel('y'),zlabel('z')
hold off

end
