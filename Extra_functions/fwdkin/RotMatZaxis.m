function [Rz] = RotMatZaxis(angle,outputsize)
%Rotates about PRINCIPAL Z-axis over angle.
%outputsize = 0 gives 3x3 matrix as output.
%outputsize = 1 gives 4x4 matrix as output.

if outputsize == 0
    
Rz = [cos(angle),  -sin(angle),  0;
      sin(angle),   cos(angle),  0;
               0,            0,  1      ];
  
end

if outputsize == 1
    
Rz = [cos(angle),  -sin(angle),  0, 0;
      sin(angle),   cos(angle),  0, 0;
               0,            0,  1, 0;
               0,            0,  0, 1   ];
    
end

end

