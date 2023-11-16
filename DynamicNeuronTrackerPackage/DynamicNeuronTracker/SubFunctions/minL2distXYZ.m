function D2 = minL2distXYZ(XI,XJ)  
% minL2distXYZ 

[m,p] = size(XJ);
px = p/3;

dtX = zeros(m, px);
dtY = zeros(m, px);
dtZ = zeros(m, px);

dtX = (XI(1:px) - XJ(:, 1:px)).^2;
dtY = (XI(px+1:2*px) - XJ(:, px+1:2*px)).^2;
dtZ = (XI(2*px+1:p) - XJ(:, 2*px+1:p)).^2;
dtXYZ = dtX + dtY + dtZ;

% min_{t}
D2 = sqrt(min(dtXYZ, [], 2));

end
