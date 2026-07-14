function [ x,y,circle ] = createCircle( radius, numberofpoints )
%Creates circle

n=numberofpoints;
th=linspace(0,2*pi,n)';
x=radius*cos(th); y=radius*sin(th);
circle=[x y]';

end

