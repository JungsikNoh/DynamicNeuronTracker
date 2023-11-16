function F = ftnGaussian3D(prmVect0, xRange, yRange, zRange)
%    Notation: x : x-position
%              y : y-position
%              z : z-position 
%              A : amplitude
%              s : x,y standard deviation
%              r : z standard deviation
%              c : background 

x0 = prmVect0(1);
y0 = prmVect0(2);
z0 = prmVect0(3);
A = prmVect0(4);
s = prmVect0(5);
r = prmVect0(6);
c = prmVect0(7);

xRange = xRange - x0;
yRange = yRange - y0;
zRange = zRange - z0;

[X,Y,Z] = meshgrid(xRange, yRange, zRange);

a = -0.5/s;
b = -0.5/r;

F = A * exp( X.^2 .* a + Y.^2 .* a + Z.^2 .*b ) + c;


end
