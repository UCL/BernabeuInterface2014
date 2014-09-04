% Uses voronoiSkel from fileexchange http://www.mathworks.co.uk/matlabcentral/fileexchange/27543-skeletonization-using-voronoi/content/voronoiSkel.m
% If the segmentation needs tweaking, there is a 'trim' parameter in voronoiSkel that can help removing spurious branches.

clear;

%%
testPlot = 0;
[fileName,dirName] = uigetfile('*.tif');

inputImg = imread(fullfile(dirName,fileName));

% Turn RGB to BW, cheat and take R channel and invert it (for Black signal
% on white background in tiff image)
plexusImg = imclose(~inputImg(:,:,1),strel('disk',5));

%%
fprintf('Getting skeleton...');
[voronoiSkeleton,vertices,edges] = voronoiSkel(plexusImg);
fprintf('finished\n');
negPlexus = ~plexusImg;
B = bwboundaries(plexusImg);
plexusBoundary = vertcat(B{:});

%%
% Produce matlab skeleton by postprocessing voronoi skeleton
skeleton = bwmorph(voronoiSkeleton,'skel','Inf');
branchpoints = bwmorph(voronoiSkeleton,'branchpoints');

%%
fprintf('Preparing plexus image...');
% [plexusR,plexusC] = find(negPlexus);
fprintf('finished\n');
fprintf('Getting distances...');
% [D,I] = pdist2([plexusC plexusR],v,'euclidean','Smallest',1);
[radiusIndex,radius] = knnsearch(plexusBoundary, vertices);
fprintf('finished\n');

[path,name,ext] = fileparts(fileName);
outputFileName = fullfile(dirName,[name '.mat']);

save(outputFileName);

if testPlot
    %% Plot plexus with circles at voronoi vertices (proportional to vessel diameter)
    radius(radius==0) = [];
    for i = 1:length(radius)
        plot(vertices(i,2),vertices(i,1),'o','MarkerSize',radius(i)/10);
        hold on;
    end
    axis ij image;
    hold off;
end

