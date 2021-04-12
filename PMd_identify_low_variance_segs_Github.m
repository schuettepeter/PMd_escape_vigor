%use this script to identify cells of low variance to exclude from
%analysis.

%first calculate variance for all cells in a folder
varianceThreshFrac = .25;

load('output_CNMF-E.mat','neuron')

for cellNum = 1:size(neuron.C,1)
    segVariance(cellNum) = var(neuron.C_raw(cellNum,:)); 
end

varThresh = max(segVariance) .* varianceThreshFrac; %overwritten below if decide to use outliers

%set this to '1' is choose to exclude outliers from variance thresh
%calculation
%remove outliers first
tempOutlier = isoutlier(segVariance,'median');
%calculate the lower bound variance threshold of cells to keep
varThresh = max(segVariance(find(~tempOutlier))) .* varianceThreshFrac;

good_neurons = segVariance > varThresh;
%sum(good_neurons) 
%find(good_neurons)

save('good_neurons.mat','good_neurons')