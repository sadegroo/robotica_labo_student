function [ xcolumn,ycolumn,zcolumn,reshapedSphere ] = createSphere( radius, numberofpoints )
%Reshapes the sphere coordinates into column vectors.

[x,y,z] = sphere(numberofpoints);
[s,d] = size(xsphere)
xcolumn = zeros;
ycolumn = zeros;
zcolumn = zeros;


for i = 1:1:s^2
    
    newx = xsphere(i);
    newy = ysphere(i);
    newz = zsphere(i);
    
    xrow(i) = newx;
    yrow(i) = newy;
    zrow(i) = newz;
    
end

xcolumn = xrow';
ycolumn = yrow';
zcolumn = zrow';
reshapedSphere = [xcolumn,ycolumn,zcolumn];
reshapedSphere = reshapedSphere*radius;

end

