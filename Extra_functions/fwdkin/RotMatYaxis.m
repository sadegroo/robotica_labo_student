function [Ry] = RotMatYaxis(angle,outputsize)
%Rotates about PRINCIPAL Y-axis over angle.
%outputsize = 0 gives 3x3 matrix as output.
%outputsize = 1 gives 4x4 matrix as output.
  
if outputsize == 0
    
Ry = [ cos(angle),   0,    sin(angle);
                0,   1,             0;
      -sin(angle),   0,    cos(angle)   ];
  
end

if outputsize == 1
    
Ry = [ cos(angle),   0,  sin(angle),  0;
                0,   1,          0,   0;
      -sin(angle),   0,  cos(angle),  0;
                0,   0,          0,   1     ];
    
end

end

