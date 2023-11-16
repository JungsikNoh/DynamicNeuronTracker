function corrfiltOut = corrfilter3D(Hscaled, B)
% corrfilter3D


%%
B = double(B);

[sH1, sH2, sH3] = size(Hscaled);
[sB1, sB2, sB3] = size(B);
n2 = numel(Hscaled(:));

Bsum0 = imfilter(B, ones(size(Hscaled)), 'full');
Bsum = Bsum0(sH1:sB1, sH2:sB2, sH3:sB3);

B2sum0 = imfilter(B.^2, ones(size(Hscaled)), 'full');
B2sum = B2sum0(sH1:sB1, sH2:sB2, sH3:sB3);
 
B_nbarxSq = n2 * (Bsum ./ n2).^2;
Bstd = sqrt((B2sum - B_nbarxSq)./(n2-1));

crsprod0 = imfilter(B, Hscaled, 'full');
crsprod = crsprod0(sH1:sB1, sH2:sB2, sH3:sB3);
corrfiltOut = crsprod ./ (n2-1) ./ Bstd; 

end
