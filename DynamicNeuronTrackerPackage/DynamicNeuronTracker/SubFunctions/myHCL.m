function [fig1, fig2, fig3, outperm, hclmap, fig4, fig5, outpermCol, hclmap2] = ...
            myHCL(inputMat, linkageMethod, varargin)
% myHCL implements Hierarchical Clustering analysis for the matrix of
% meanActivityOfNeuronsPerStimulus.


%% input

ip = inputParser;
ip.addParameter('rowNames', '')
ip.addParameter('xlabel', 'X')
ip.addParameter('ylabel', 'Y')
ip.addParameter('figFlag', 'on')
ip.addParameter('title', '')

ip.parse(varargin{:});
p = ip.Results;

%% input heatmap

fig1 = figure('Visible', p.figFlag);
q = quantile(inputMat(:), [0.005, 0.995]);
imagesc(inputMat, [q(1), q(2)])
colorbar; colormap(jet)

xlabel(p.xlabel)
ylabel(p.ylabel)
title(p.title)

if ~isempty(p.rowNames)
    ax = gca;
    ax.YTick = 1:size(inputMat, 1);
    ax.YTickLabel = p.rowNames;
    ax.YTickLabelRotation = 0;
end
 


%% hcl

D = pdist(inputMat);
tree = linkage(inputMat, linkageMethod);
leafOrder = optimalleaforder(tree, D);

fig2 = figure('Visible', p.figFlag);
[H,T,outperm] = dendrogram(tree, size(inputMat, 1), 'Reorder', leafOrder, ...
    'Orientation', 'left', 'Labels', p.rowNames);

xlabel('Distance')
ylabel(p.ylabel)

%% hclmap

hclmap = inputMat(outperm, :);

fig3 = figure('Visible', p.figFlag);
imagesc(hclmap, [q(1), q(2)])
colorbar; colormap(jet)

xlabel(p.xlabel)
ylabel(p.ylabel)
title(p.title)

if ~isempty(p.rowNames)
    ax = gca;
    ax.YTick = 1:size(hclmap, 1);
    ax.YTickLabel = p.rowNames(outperm);
    ax.YTickLabelRotation = 0;
end

%% hcl columns for better visualization
 
inputMat = hclmap';
D = pdist(inputMat);
tree = linkage(inputMat, linkageMethod);
leafOrder = optimalleaforder(tree, D);

fig4 = figure('Visible', p.figFlag);
[H,T,outpermCol] = dendrogram(tree, size(inputMat, 1), 'Reorder', leafOrder, ...
    'Orientation', 'left');

xlabel('Distance')
ylabel(p.xlabel)

%% hclmap2

hclmap2 = hclmap(:, outpermCol);

fig5 = figure('Visible', p.figFlag);
imagesc(hclmap2, [q(1), q(2)])
colorbar; colormap(jet)

xlabel(p.xlabel)
ylabel(p.ylabel)
title(p.title)

if ~isempty(p.rowNames)
    ax = gca;
    ax.YTick = 1:size(hclmap, 1);
    ax.YTickLabel = p.rowNames(outperm);
    ax.YTickLabelRotation = 0;
end

end
