function D2 = nanhamdistXYZ_ignoreAllnan(XI,XJ)  
% nanhamdistXYZ_ignoreAllnan computes Hamming distance ignoring coordinates with NaNs

[m,p] = size(XJ);
px = p/3;
nesum = zeros(m,1);
pstar = zeros(m,1);
for q = 1:px
    notnanx = ~(isnan(XI(q)) & isnan(XJ(:,q)));
    notnany = ~(isnan(XI(q+px)) & isnan(XJ(:,q+px)));
    notnanz = ~(isnan(XI(q+2*px)) & isnan(XJ(:,q+2*px)));
    nex = ((XI(q) ~= XJ(:,q)) & notnanx);
    ney = ((XI(q+px) ~= XJ(:,q+px)) & notnany);
    nez = ((XI(q+2*px) ~= XJ(:,q+2*px)) & notnanz);
    ne = max(max(nex, ney), nez);
    notnan = notnanx .* notnany .* notnanz;
    nesum = nesum + ne;
    pstar = pstar + notnan;
end
D2 = nesum./pstar; 
