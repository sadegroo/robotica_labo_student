function [ CORD ] = ConvertSphere( x,y,z,JacobianColumns )
%takes the sphere coordinates and arranges them in a 3 row matrix,
%      [
%        x
%        y
%        z
%          ]

[a,b] = size(x);
[c,d] = size(y);
[e,f] = size(z);

xsize = a*b;
ysize = c*d;
zsize = e*f;

X = reshape(x,[1,xsize]);
Y = reshape(y,[1,ysize]);
Z = reshape(z,[1,zsize]);
dummy = zeros(1,xsize);

if JacobianColumns == 1
CORD = [X];
end

if JacobianColumns == 2
CORD = [X;Y];
end

if JacobianColumns == 3
CORD = [X;Y;Z];
end

if JacobianColumns == 4
CORD = [X;Y;Z;dummy];
end

if JacobianColumns == 5
CORD = [X;Y;Z;dummy;dummy];
end

if JacobianColumns == 6
CORD = [X;Y;Z;dummy;dummy;dummy];
end

end

