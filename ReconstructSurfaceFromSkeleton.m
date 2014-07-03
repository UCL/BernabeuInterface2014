% You will have to make sure that MATLAB can find the VMTK dynamic libraries at runtime.

function [] = ReconstructSurfaceFromSkeleton(filename, pixelsPerUm)

    % Depending on whether we assume that the vessels remain circular in cross
    % section (1.0) or they collapse after fixation (2/pi).
    radiiFudgeFactor = 1.0;
    %radiiFudgeFactor = 2 / pi;

    % Extract datasetName to be used as identifier
    [~,datasetName,ext] = fileparts(filename);
    assert(strcmp(ext, '.mat'), 'Wrong file extension, it should be ''.mat''. Use SkeletonizeTiffPlexus to process a ''.tif'' before calling the current function.')
    
    % Load dataset
    load(filename);
    
    pixelToUm = 1.0 / pixelsPerUm;

    % Choose according to the variable containing the radii information
    %radiiVariable = radii;
    radiiVariable = radius;

    % Run writer
    lengthRadiusForEachEdge = VTP_writer(vertices, edges, radiiVariable, datasetName, pixelToUm, radiiFudgeFactor);

    % Histogram of number of segments for a given diameter
    hist(2*radiiFudgeFactor*radiiVariable*pixelToUm, 100)
    xlabel('Diameter (um)')
    ylabel('Number of network segments')
    %title('Network diameter histogram')
    print('-dpng', [datasetName '.png'])

    % Histogram of length covered by segments of a given diameter
    largestRadius = ceil(max(lengthRadiusForEachEdge(:,2)))
    range = 0+eps:0.05:largestRadius;
    lengths(1:length(range)) = 0;
    for edgeId=1:size(lengthRadiusForEachEdge,1),
        for i=1:length(range)-1,
            if (lengthRadiusForEachEdge(edgeId,2) >= range(i) && lengthRadiusForEachEdge(edgeId,2) < range(i+1)),
                lengths(i) = lengths(i) + lengthRadiusForEachEdge(edgeId,1);
            end
        end
    end
    
    % Do some stats
    parmhat = lognfit(2*range, [], [], lengths);    
    mu = parmhat(1);
    sigma = parmhat(2);
    mean = exp(mu + sigma^2/2);
    variance = exp(2*mu + sigma^2)*(exp(sigma^2) - 1);
    mode = exp(mu - sigma^2);
    median = exp(mu);
    
    
    
    figure; bar(2*range,lengths)
    axis([0 40 0 800])
    set(gca,'XTick', 0:2:40)
    xlabel('Diameter (um)')
    ylabel('Total length covered (um)')
    %title('Network radii histogram')

    % Plot a vertical line to mark the mode of the lognormal distribution
    hold on 
    line = plot([mode mode], ylim, 'r')
    legend_txt= sprintf('mode = %.2f um', mode);
    legend(line, legend_txt);

    set(gca, 'FontSize', 17)
    set(findall(gcf, 'type', 'text'), 'FontSize', 17)
    
    print('-dpng', [datasetName '_length.png'])
    
    total_network_length = sum(lengths)
        
    % Compute the value of delta x required to satisfy that no more
    % than fractionToBeCovered of the total length of the network has fewer
    % than diameterInLatticeUnits lattice sites accross.
    diameterInLatticeUnits = 7;
    totalLength = sum(lengths);
    partialSum = 0;
    fractionToBeCovered = 0.05;
    for i=1:length(lengths),
        partialSum = partialSum + lengths(i);
        if partialSum > (fractionToBeCovered*totalLength),
            sprintf('%f length covered at %f um', fractionToBeCovered, range(i))
            delta_x = 2*range(i)/(diameterInLatticeUnits-1);
            sprintf('Set delta x to %f um', delta_x)
            break        
        end
    end

    % based on the fact that we want nu_inf to map tau=0.6
    nu_inf = 3.85e-6;
    delta_t = (delta_x*1e-6)^2 * 0.1 / (3 * nu_inf);
    sprintf('Set delta t to %e s', delta_t)

    networkSize = (max(vertices) - min(vertices)) * pixelToUm * 1e-6;
    longestPath = 2*networkSize(1) + networkSize(1)/2;
    t_diff = (2 * largestRadius * 1e-6)^2 / nu_inf
    sprintf('Run %f timesteps', t_diff/delta_t)

    % Using vmtkcenterlinemodeller to generate the surface from the skeleton takes forever. I checked with Luca 
    % and vmtkpolyballmodeller is a valid option as long as the network nodes are quite close together (which is our case)
    vmtkPipeline = ['vmtkpolyballmodeller -ifile ' datasetName '.vtp -radiusarray Radius -dimensions 1024 1024 256 --pipe vmtkmarchingcubes -ofile ' datasetName '_tubed.stl']
    system(vmtkPipeline);

    % Filter configuration values 0.1 and 30 taken from VMTK tutorial "PREPARE A SURFACE FOR MESH GENERATION"
    vmtkPipeline = ['vmtksurfacesmoothing -ifile ' datasetName '_tubed.stl' ' -passband 0.1 -iterations 30 -ofile ' datasetName '_tubed_smoothed.stl']
    system(vmtkPipeline);

    close all;

end
