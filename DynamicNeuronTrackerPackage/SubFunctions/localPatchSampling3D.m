function [avgPatch, Patches] = localPatchSampling3D(x,y,z,f, ...
        bandWidthX, bandWidthZ, imgArray)
% localPatchSampling3D

%%

winLengthX = 2*bandWidthX + 1;
winLengthZ = 2*bandWidthZ + 1;
Patches = nan(winLengthX, winLengthX, winLengthZ, numel(f));

for i=1:numel(f)    
        xwin = x(i)-bandWidthX:x(i)+bandWidthX;
        ywin = y(i)-bandWidthX:y(i)+bandWidthX;
        zwin = z(i)-bandWidthZ:z(i)+bandWidthZ;
        patchtmp = imgArray(ywin, xwin, zwin, f(i));
        Patches(:, :, :, i) = patchtmp;
end
    
avgPatch = mean(Patches, 4, 'omitnan');

end

