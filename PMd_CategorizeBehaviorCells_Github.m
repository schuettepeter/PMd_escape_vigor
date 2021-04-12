%% The following code was used to determine which putative neurons significantly 
%encoded defensive behaviors for the manuscript "Dorsal premammillary hypothalamic 
%projection to periaqueductal gray controls escape vigor from innate and
%conditioned threats."

%Briefly, the neural data is fit to the behavioral data by Generalized
%Linear model. For each recorded putative cell, the actual coefficient is
%compared to a distribution of coefficients calculated by rolling the
%neural data in relation to the behavioral data.  If the actual coefficient
%is greater than or less that 95% of this permuted distribution, then the
%given cell significantly encodes the behavior in question.

clearvars -except folders mouseNum assayNum

numSamplesPre = round(7.5 .* 3); %fps * # seconds

%load a matrix of neural data ('neuron')
load('output_CNMF-E.mat', 'neuron')
load('Tracking.mat'); 
load('good_neurons.mat')

%load the relevant behavioral data. (vectors of 0's and 1's, the same
%length as the neural data)
if assayNum==1 | assayNum==2 | assayNum==3 | assayNum==4
load('BehaviorMS.mat','approachIndicesMS','approachFrameMS', 'stretchIndicesMS','stretchFrameMS','escapeFrameMS','escapeIndicesMS','freezeFrameMS','freezeIndicesMS')  
numBehavs = 4;
behavs = {'approach','stretch','escape','freeze'};
end
good_neurons = find(good_neurons);
%%

for behavNum = 1:numBehavs
    
            if assayNum==1 | assayNum==2 | assayNum==3 | assayNum==4
                if behavNum==1
                        behavFrame = approachFrameMS;
                elseif behavNum==2
                        if exist('stretchFrameMS')
                            behavFrame = stretchFrameMS;
                        else
                            continue
                        end
                elseif behavNum==3
                        behavFrame = escapeFrameMS;
                elseif behavNum==4
                        if exist('freezeFrameMS')
                            behavFrame = freezeFrameMS;
                        else
                            continue
                        end
                end
            end
            
            behavIndices = zeros(length(neuron.C_raw),1);  
                                        
                behavFrame(:,2) = behavFrame(:,2);
                behavFrame(:,1) = behavFrame(:,1) - numSamplesPre;
                
                behavFrame(find(behavFrame<1)) = 1;
                behavFrame(find(behavFrame>length(behavIndices))) = length(behavIndices);
                                
            for behavCount = 1:size(behavFrame,1)
                behavIndices([behavFrame(behavCount,1):behavFrame(behavCount,2)]) = 1;
            end
                        
            Indices = behavIndices;
                     
        for seg = 1:length(good_neurons)

            CalciumTrace = neuron.C_raw(good_neurons(seg),:)';
            
            while length(Indices) > length(CalciumTrace)
                Indices = Indices(1:end-1);
            end
            while length(Indices) < length(CalciumTrace)
                Indices = [Indices; 0];
            end
            
            tbl = table(Indices, CalciumTrace);
            lm = fitglm(tbl, 'CalciumTrace~Indices');

            coeffSeg(good_neurons(seg),:) = table2array(lm.Coefficients(2:end,[1,4]));

        end
    
%create bootstrap distributions for each neuron to compare with output from previous section -- 
% ROLL THE CALCIUM DATA, RATHER THAN SHUFFLE BEHAVIOR EPOCHS

    %determine the number of samples to use (# of shuffles)
    samples = floor(length(CalciumTrace) ./ 45) - 10;

    disp(['building bootstrap distribution'])
    
    InOut = behavFrame;
    
    for iter = 1:samples

            CalciumTrace = neuron.C_raw(seg,:)';
            CalciumTrace = [CalciumTrace; CalciumTrace(1:round(45 .* iter)-1)]; %add a bit to the end
            CalciumTrace = CalciumTrace(round(45 .* iter):end); %remove a bit from the beginning        
        
            for seg = 1:length(good_neurons)
                tbl = table(Indices, CalciumTrace);
                lm = fitglm(tbl, 'CalciumTrace~Indices');
                coeffSeg_shuff(good_neurons(seg),iter) = table2array(lm.Coefficients(2:end,1));
            end
    end
    
%determine, for each seg, where it's coefficient falls in the distribution
%-- EITHER ABOVE OR BELOW (so either positively or negatively modulated by behavior)
    for seg = 1:length(good_neurons)
        %positively modulated by behavior
        temp = length(find(coeffSeg_shuff(good_neurons(seg),:) > coeffSeg(good_neurons(seg),1))) ./ samples;
        if temp <= .05 & coeffSeg(good_neurons(seg),2) <= .05 %p-value of bootstrap and coefficient.
           behavSeg(good_neurons(seg)) = 1;
        else
           behavSeg(good_neurons(seg)) = 0;
        end
        
        %negatively modulated by behavior
        temp = length(find(coeffSeg_shuff(good_neurons(seg),:) < coeffSeg(good_neurons(seg),1))) ./ samples; %fraction bootstrapped vals that are LESS THAN ACTUAL VALUE
        if temp <= .05 & coeffSeg(good_neurons(seg),2) <= .05 %p-value of bootstrap and coefficient.
           behavSeg(good_neurons(seg)) = -1;
        end
    end 

    behavSegAll{behavNum} = behavSeg;
    coeffSegAll{behavNum} = coeffSeg;
    coeffSeg_shuffAll{behavNum} = coeffSeg_shuff;
    
    clearvars behavSeg coeffSeg coeffSeg_shuff
end

%Save results with each recording session's data.
save('Seg.mat','behavSegAll','coeffSegAll','coeffSeg_shuffAll','behavs')       