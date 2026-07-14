%% Voorbeeld 4.1 uit Modern Robotics (Lynch & Park): 3R ruimtelijke open keten
% Voorwaartse kinematica van de 3R-keten uit figuur 4.3, symbolisch berekend
% met de product-of-exponentials formule in ruimteframe-vorm:
%
%   T(theta) = e^[S1]th1 * e^[S2]th2 * e^[S3]th3 * M
%
% De berekening gebruikt de Modern Robotics MATLAB-bibliotheek
% (modernrobotics/packages/MATLAB/mr) samen met de Symbolic Math Toolbox.
%% Symbolische variabelen
clear; clc;
syms L1 L2 positive
syms th1 th2 th3 real

%% Gegevens uit voorbeeld 4.1 (figuur 4.3)
% Home-configuratie van het eindeffectorframe {3} t.o.v. vast frame {0}
M = [ 0  0  1  L1;
      0  1  0  0;
     -1  0  0 -L2;
      0  0  0  1];

% Schroefassen S_i = (omega_i, v_i) in het vaste (ruimte)frame:
%   i |  omega_i   |  v_i
%   1 | (0, 0, 1)  | (0,  0,   0)
%   2 | (0,-1, 0)  | (0,  0, -L1)
%   3 | (1, 0, 0)  | (0,-L2,   0)
S1 = [0;  0; 1; 0;   0;   0];
S2 = [0; -1; 0; 0;   0; -L1];
S3 = [1;  0; 0; 0; -L2;   0];
Slist = [S1, S2, S3];

thetalist = [th1; th2; th3];

%% Voorwaartse kinematica met de MR-bibliotheek
T = FKinSpace(M, Slist, thetalist);
T = simplify(T);
% T = simplify(T, 'Steps', 100) % gebruik evt meer Steps in simplify
% T = collect(T, [symexpr1, symexpr2, ...]); % om termen te % factoriseren

disp('T(theta) volgens FKinSpace:')
disp(T)

%% Omvormen naar matlab functie
% 1) mechanische parameters invullen
L1_val = 1;  L2_val = 1;
T_sub = subs(T, [L1 L2], [L1_val L2_val]);

% 2) omzetten naar function handle met één vectorargument
%    (de accolades rond de vector zorgen dat theta als 1 vector binnenkomt)
fkin = matlabFunction(T_sub, 'Vars', {thetalist});

% gebruik:
T_num = fkin([0; pi/2; pi])

%% Plot frame
plotFrame(eye(3)*0.6,zeros(3,1))
hold on
plotFrame(T_num,0)
hold off
view([46.667 42.667])
axis equal

% output naar latex (; wissen)
mat2latex(T, shorthand=true, math="display");
mat2latex(round(T_num,2), shorthand=true, math="display");
