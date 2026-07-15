function T = DH(alpha,a,d,theta)
% Compute a single DH-transformation: takes a single row of DH-table and computes transformation matrices
% According to the classic Denavit-Hartenberg method.
% Arguments:
%   one row of DH table in order: alpha - a - d - theta
%   alpha and theta should be entered in RADIANS!
% Return Values:
%   T: Returns a 4x4 double DH-transform

% Generate DH matrix
T =         [cos(theta)  -sin(theta)*cos(alpha)   sin(theta)*sin(alpha)  a*cos(theta) ;
             sin(theta)   cos(theta)*cos(alpha)  -cos(theta)*sin(alpha)  a*sin(theta) ;
             0            sin(alpha)              cos(alpha)             d            ;
             0            0                       0                      1            ];
            
end