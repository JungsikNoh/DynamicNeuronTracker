function outMIPmat = simpleVol2MIP(vol)
% simpleVol2MIP convert a small volume image into a MIP 2D-image.


[s1, s2, s3] = size(vol);
outMIPmat = nan(s1+1+s3, s2+1+s3);

pjz = max(vol, [], 3);
pjy0 = squeeze(max(vol, [], 1));
pjx = squeeze(max(vol, [], 2));
pjy = transpose(pjy0);

outMIPmat(1:s1, 1:s2) = pjz;
outMIPmat(s1+1+(1:s3), 1:s2) = pjy;
outMIPmat(1:s1, s2+1+(1:s3)) = pjx;

end
