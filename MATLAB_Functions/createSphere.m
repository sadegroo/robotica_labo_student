function [ xcolumn,ycolumn,zcolumn,reshapedSphere ] = createSphere( radius,numberofpoints )
%Creates sphere with given number of points and radius.
%Reshapes the sphere into a column vector with 3 columns (x,y,z
%coordinates).

[x,y,z] = sphere(numberofpoints);
[s,d] = size(x);

for i = 1:1:s^2
    
    newx = x(i);
    newy = y(i);
    newz = z(i);
    
    xrow(i) = newx;
    yrow(i) = newy;
    zrow(i) = newz;
    
end

xcolumn = radius*xrow';
ycolumn = radius*yrow';
zcolumn = radius*zrow';
reshapedSphere = [xcolumn,ycolumn,zcolumn];
%reshapedSphere = reshapedSphere*radius;

end



