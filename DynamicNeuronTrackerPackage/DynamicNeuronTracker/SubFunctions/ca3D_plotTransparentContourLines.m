function ca3D_plotTransparentContourLines(C, alpha, colvec)
% ca3D_plotTransparentContourLines
% Transparent contour lines of ROIs greatly help visualize ROIs on top of a
% MIP video.


numPts1 = C(2, 1);
pts1 = C(:, 2:1+numPts1);
numPts2 = C(2, 2+numPts1);
pts2 = C(:, 3+numPts1:2+numPts1+numPts2);
numPts3 = C(2, 3+numPts1+numPts2);
pts3 = C(:, 4+numPts1+numPts2:3+numPts1+numPts2+numPts3);

l1 = plot(pts1(1,:), pts1(2, :), 'Color', [colvec, alpha], 'LineWidth', 1);
l2 = line(pts2(1,:), pts2(2, :), 'Color', [colvec, alpha], 'LineWidth', 1);
l3 = line(pts3(1,:), pts3(2, :), 'Color', [colvec, alpha], 'LineWidth', 1);

end