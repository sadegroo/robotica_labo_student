function [ Tfull,Tparts,Tcumul ] = DH_modified( DH_table )
% Compute DH-transformations: takes every row of DH-table and computes transformation matrices
% According to the MODIFIED Denavit-Hartenberg method (Craig, "Introduction
% to Robotics: Mechanics and Control", eq. 3.6).
% Arguments:
%   DH_table: n x 4 double or sym with columns in the order:
%       alpha_(i-1) - a_(i-1) - d_i - theta_i
%       alpha and theta should be entered in RADIANS!
% Return Values:
%   Tfull: Returns a 4x4 double or sym DH-transform of the full transform from frame 0 to the
%   last frame
%   Tparts: Returns a cell array with the 4x4 double DH transforms for each row of the DH table
%       cell i contains T(i-1)->i. To address matrix, type >> Tparts{i}
%   Tcumul: Returns a cell array of the intermediate transforms from frame 0
%       to the frame corresponding to the respective row
%       cell i contains T0->i. To address matrix, type >> Tcumul{i}
%
% MODIFIED DH CONVENTION (differs from the classical convention of DH_full):
%   ^{i-1}T_i = Rot(x, alpha_(i-1)) * Trans(x, a_(i-1))
%             * Rot(z, theta_i)     * Trans(z, d_i)
%
%             | c_t       -s_t        0      a     |
%             | s_t*c_a    c_t*c_a   -s_a   -d*s_a |
%             | s_t*s_a    c_t*s_a    c_a    d*c_a |
%             |  0          0         0      1     |
%
%   i.e. the x-axis screw uses the PREVIOUS link's alpha and a, and is
%   applied BEFORE the z-axis screw (classical DH applies it after).

% initialize
temp = eye(4);
sz = size(DH_table,1);
Tparts = cell(1,sz);
Tcumul = cell(1,sz);

for i = 1:sz
% address a line of the DH table
alpha = DH_table(i,1);
a = DH_table(i,2);
d = DH_table(i,3);
theta = DH_table(i,4);

% Generate modified DH matrix (Craig eq. 3.6)
Tparts{i} = [cos(theta)             -sin(theta)             0            a             ;
             sin(theta)*cos(alpha)   cos(theta)*cos(alpha) -sin(alpha)  -d*sin(alpha)  ;
             sin(theta)*sin(alpha)   cos(theta)*sin(alpha)  cos(alpha)   d*cos(alpha)  ;
             0                       0                      0            1             ];

%right multiply with new DH matrix
temp = temp*Tparts{i};

% Return intermediate transform up to frame "i"
Tcumul{i} = temp;
end

% Return full DH matrix
Tfull = temp;

end
