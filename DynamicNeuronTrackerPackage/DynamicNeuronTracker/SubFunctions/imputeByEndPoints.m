function impVec = imputeByEndPoints(rowvec)
% imputeByEndPoints

x = rowvec;
x(isnan(x)) = 0;

tmp = rle(x);
l0 = numel(tmp);

vals = tmp(1:2:l0-1);
lens = tmp(2:2:l0);

vals1 = vals;
if (tmp(1) == 0)
    vals1 = [vals1(2), vals1];
end

ind = find(vals1 == 0);
ind1 = ind - 1;
endingValues = vals1(ind1);

impVals = vals1;
impVals(ind) = endingValues;

if (tmp(1) == 0)
    impVals = impVals(2:end);
end

mat0 = [impVals; lens];
vec0 = mat0(:)';

impVec = irle(vec0);

end
