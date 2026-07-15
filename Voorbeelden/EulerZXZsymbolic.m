%% Voorbeeld ZXZ Eulerhoek-parametrisatie
clear; clc;
syms psi theta varphi

% Antisymmetrische matrixvoorstellingen van rotatievectoren
omZ1 = VecToso3([0;0;psi])
omX2 = VecToso3([theta;0;0])
omZ3 = VecToso3([0;0;varphi])

% individuele rotatiematrices
RotZ1 = MatrixExp3(omZ1)
RotX2 = MatrixExp3(omX2)
RotZ3 = MatrixExp3(omZ3)

% Product van matrixexponenten (in juiste volgorde!)
rotZXZ = RotZ1*RotX2*RotZ3

% of in 1 keer met deze functie (>>help Euler2RotMat)
rotZXZ_alt = Euler2RotMat([psi theta varphi],'zxze')

% Beide manieren moeten hetzelfde resultaat opleveren
klopt_het = isequal(rotZXZ,rotZXZ_alt) % 1=true, 0=false

% Omvormen naar numerieke functie handle (later handig bij voorwaarste kinematica)
f_rotZXZ = matlabFunction(rotZXZ, 'Vars', {[psi theta varphi]});

% Evalueer de numerieke functie met hoeken [1, 2, 1]
rotZXZ_num = f_rotZXZ([1,2,1])

% Inverse kinematisch algoritme om Eulerhoeken te recupereren. 
% (>>help RotMat2Euler)
[eulzxz1,eulzxz2] = RotMat2Euler(rotZXZ_num,'zxze')

% 2 antwoorden?! eulzxz1 is wat we ingaven, wat met eulzxz2?
rotZXZ_num2 = Euler2RotMat(eulzxz2,'zxze');

verschil = rotZXZ_num - rotZXZ_num2 % zéér klein, stelt dezelfde rotatie voor!
norm_vershil = norm(verschil) % Frobenius norm = matrix norm (wortel van som van alle abs(a_ij)^2)

% we kunnen eender welke Eulerparametrisatie opvragen, dit geeft natuurlijk
% andere hoeken
[eulzyx1,eulzyx2] = RotMat2Euler(rotZXZ_num,'zyxe')

% Of de as-hoek parametrisatie (>>help RotMat2AxisAngle)
[K, theta] = RotMat2AxisAngle(rotZXZ_num)


%% Print to latex
mat2latex(rotZXZ, shorthand=true, math="display")