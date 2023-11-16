function [Xoutlt, Youtlt, Zoutlt, meanIntVec, indVec] = ...
    selectTracks_Unique_WithMinNumBrightPS(Xmat, Ymat, Zmat, brPSMat, imgArray, minNumBrightPS)
% selectTracks_Unique_WithMinNumBrightPS


disp('== num of input tracks')
disp(size(Xmat, 1))
indVec0 = reshape(1:size(Xmat, 1), [], 1);

%% make output tracks unique

tracksXYZ1_out = [Xmat, Ymat, Zmat];
 
tracksXYZ1_tmp = tracksXYZ1_out;
tracksXYZ1_tmp(isnan(tracksXYZ1_out)) = 0;

[~, ia, ~] = unique(tracksXYZ1_tmp, 'rows', 'stable');
%size(C)
Xmatu = Xmat(ia, :);
Ymatu = Ymat(ia, :);
Zmatu = Zmat(ia, :);
indVec1 = indVec0(ia);

% check an all-nan track
indic = (all(isnan(Xmatu), 2));
disp('== index of the all-nan track')
disp(find(indic))
Xmatu = Xmatu(~indic, :);
Ymatu = Ymatu(~indic, :);
Zmatu = Zmatu(~indic, :);
indVec2 = indVec1(~indic);
disp('== num of unique tracks')
disp(size(Xmatu, 1))

%% tracks consisting of only bright point sources (in 1 pixel distance)

Xmatu_fir = nan(size(Xmatu));
Ymatu_fir = nan(size(Ymatu));
Zmatu_fir = nan(size(Zmatu));

for itr = 1:size(Xmatu, 1)
    
    xvec = Xmatu(itr, :);
    yvec = Ymatu(itr, :);
    zvec = Zmatu(itr, :);
    ptInd = nan(1, numel(xvec));

    for fr = 1:numel(xvec)
        if ~isnan(yvec(fr)) 
            psArea = brPSMat(yvec(fr)-1:yvec(fr)+1, xvec(fr)-1:xvec(fr)+1,zvec(fr)-1:zvec(fr)+1, fr);
            if (sum(psArea(:)) == 1)            % select clear firing events
                ptInd(fr) = 1;
            else
                ptInd(fr) = sum(psArea(:));
            end
        end
    end
    
    indIntersect = (ptInd == 1);            % clear/Bright firing events
    x = xvec; y=yvec; z=zvec;
    x(~indIntersect) = NaN;
    y(~indIntersect) = NaN;
    z(~indIntersect) = NaN;
    Xmatu_fir(itr,:) = x;
    Ymatu_fir(itr,:) = y;
    Zmatu_fir(itr,:) = z;
end

%% make sure the uniqueness

tracksXYZ1_out = [Xmatu_fir, Ymatu_fir, Zmatu_fir];
 
tracksXYZ1_tmp = tracksXYZ1_out;
tracksXYZ1_tmp(isnan(tracksXYZ1_out)) = 0;

[~, ia, ~] = unique(tracksXYZ1_tmp, 'rows', 'stable');
%size(C)
Xmatu_fir2 = Xmatu_fir(ia, :);
Ymatu_fir2 = Ymatu_fir(ia, :);
Zmatu_fir2 = Zmatu_fir(ia, :);
indVec3 = indVec2(ia);

% check an all-nan track
indic = (all(isnan(Xmatu_fir2), 2));
disp('== index of the all-nan track')
disp(find(indic))
Xout = Xmatu_fir2(~indic, :);
Yout = Ymatu_fir2(~indic, :);
Zout = Zmatu_fir2(~indic, :);
indVec4 = indVec3(~indic);
disp('== num of unique brPS tracks')
disp(size(Xout, 1))
%indVec = indVec4;

%% tracks with minimum number of bright point sources

lt0 = sum(~isnan(Xout), 2);
disp('== summary of lengths of the period when bright PSs satisfy corrThreshold')
disp(summaryStatistics(lt0))
%disp(quantile(lt0, (1:10)./10))

% exclude too short ones
indlt = (lt0 >= minNumBrightPS);
ratio0 = sum(indlt)/numel(indlt);
disp('== proportion of the tracks satisfying minFiringFramesOfNeurons')
disp(ratio0)

Xoutlt = Xout(indlt, :);
Youtlt = Yout(indlt, :);
Zoutlt = Zout(indlt, :);
indVec5 = indVec4(indlt);

indVec = indVec5;

%%  double-check brightness 

maxIntVec = nan(size(Xoutlt, 1), 1);
meanIntVec = nan(size(Xoutlt, 1), 1);

for k=1:size(Xoutlt, 1)
    x = round(Xoutlt(k,:));
    y = round(Youtlt(k,:));
    z = round(Zoutlt(k,:));
    %f = find(~isnan(Xmatu(k,:)));
    tmp = nan(1, numel(x));
    for i=1:numel(x)
        if ~isnan(x(i))
            tmp(i) = imgArray(y(i),x(i),z(i),i);
        end
    end
    %tracks_Bright3(k).maxInt = max(tmp(:), [], 'omitnan');
    maxIntVec(k) = max(tmp(:), [], 'omitnan');
    meanIntVec(k) = nanmean(tmp(:));
end

%% report

disp('== num of unique tracks with minimum number of bright point sources that satisfy corrThreshold')
disp(size(Xoutlt, 1))

end
