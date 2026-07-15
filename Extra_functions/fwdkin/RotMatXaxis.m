function [Rx] = RotMatXaxis(angle,outputsize)
%Rotates about PRINCIPAL X-axis over angle.
%outputsize = 0 gives 3x3 matrix as output.
%outputsize = 1 gives 4x4 matrix as output.

if outputsize == 0
    
Rx = [1,            0,           0;
      0,   cos(angle), -sin(angle);
      0,   sin(angle),  cos(angle)      ];
  
end

if outputsize == 1
    
Rx = [1,            0,           0, 0;
      0,   cos(angle), -sin(angle), 0;
      0,   sin(angle),  cos(angle), 0;
      0,            0,           0, 1   ];
    
end

end

