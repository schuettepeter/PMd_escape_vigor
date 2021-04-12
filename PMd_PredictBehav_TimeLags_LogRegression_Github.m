%% THIS IS THE SCRIPT USED TO PREDICT BEHAVIORS, USING 5-FOLD LOGISTIC REGRESSION, 
%AS DESCRIBED IN THE MANUSCRIPT "Dorsal premammillary hypothalamic projection to 
%periaqueductal gray controls escape vigor from innate and conditioned threats"

clear all

%SET THE FILEPATH FOR EACH RECORDED SESSION

%Toy rat
    folders{1,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\611\1';
    folders{2,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\612\1';
    folders{3,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\613\1';
    folders{4,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\616\1';
    folders{5,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\618\1';
    folders{6,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\620\1';
    folders{7,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\621\1';
    folders{8,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\665\1';
    folders{9,1} = 'G:\PMD_Miniscope_Rat_Shock\ToyRat\668\1';
    
%Simple rat
    folders{1,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\611\1';
    folders{2,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\612\1';
    folders{3,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\613\1';
    folders{4,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\616\1';
    folders{5,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\618\1';
    folders{6,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\620\1';
    folders{7,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\621\1';
    folders{8,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\665\1';
    folders{9,2} = 'G:\PMD_Miniscope_Rat_Shock\Rat\668\1';
    
%Shock Habituation
    folders{1,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\611\1';
    folders{2,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\612\1';
    folders{3,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\613\1';
    folders{4,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\616\1';
    folders{5,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\618\1';
    folders{6,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\620\1';
    folders{7,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\621\1';
    folders{8,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\665\1';
    folders{9,3} = 'G:\PMD_Miniscope_Rat_Shock\Preshock\668\1';
    
%Fear retrieval
    folders{1,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\611\1';
    folders{2,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\612\1';
    folders{3,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\613\1';
    folders{4,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\616\1';
    folders{5,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\618\1';
    folders{6,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\620\1';
    folders{7,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\621\1';
    folders{8,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\665\1';
    folders{9,4} = 'G:\PMD_Miniscope_Rat_Shock\Postshock\668\1';       
   
%% load necessary data into struct

%First collect all relevant data into a struct including vectors of
%behavioral instances and neural data
for mouseNum = 1:size(folders,1)
    for assayNum = 1:size(folders,2)
   
        if isempty(folders{mouseNum,assayNum})
            continue
        end
        
        cd(folders{mouseNum,assayNum})
        
        load('output_CNMF-E.mat','neuron'); load('good_neurons.mat')
        nrn_all = neuron.C_raw;
        nrn = neuron.C_raw(find(good_neurons),:);
        
        load('BehaviorMS.mat','approachIndicesMS','stretchIndicesMS','escapeIndicesMS','freezeIndicesMS');


           behavs = {'approachIndicesMS','stretchIndicesMS','escapeIndicesMS','freezeIndicesMS'};
           behavIndices = [approachIndicesMS';stretchIndicesMS';escapeIndicesMS';freezeIndicesMS'];            
        
        %adjust if diff length than neural data
        while length(behavIndices) < length(nrn)
            behavIndices = [behavIndices, behavIndices(:,end)];
        end
        while length(behavIndices) > length(nrn)
            behavIndices = behavIndices(:,1:end-1);
        end
        
        DataAll{1}{mouseNum,assayNum} = nrn_all;        
        DataAll{2}{mouseNum,assayNum} = nrn;
        DataAll{3}{mouseNum,assayNum} = behavIndices;
        DataAll{4}{mouseNum,assayNum} = behavs;
        
        clearvars nrn closedArmIndicesMS openArmIndicesMS headDipIndicesMS approachIndicesMS stretchIndicesMS escapeIndicesMS freezeIndicesMS
        
    end
end

DataAllFieldIds = {'all neural activity','good neural activity','behavior indices','behavior order', 'EPM to Rat activity', 'EPM to Toy Rat activity','Toy Rat to Rat activity'};

%% Is neural activity predictive of behaviors at behavior onset, as well as in the future?  

usePCA = 1;
minVarAcc = 80;
timeLags = round([0 2 4 6 8 10] .* 7.5); %amount of seconds to shift training data pre-escape.
sampleDur = 2; % length of training and testing data for each behavior in seconds
sampleDur = round(sampleDur .* 7.5); %convert 'sampleDur' to frames (at 7.5fps).
assayNum = 4; %Set the assay # here-- 1:toy 2:rat 3:shock hab 4:shock ext
behavNumber = 2; %Set the behavior # here--1:approach, 2:stretch, 3: escape, 4: freeze

    for mouseNum = 1:size(folders,1)
        cd(folders{mouseNum,assayNum})
        load('good_neurons.mat')
        escapeCellPosId{mouseNum} = find(good_neurons == 1); 
        escapeCellNegId{mouseNum} = find(good_neurons == 1); 
    end

%create empty cell matrix for each mouse x each time lag.
escapePosDff = cell(size(folders,1),length(timeLags));
escapeNegDff = cell(size(folders,1),length(timeLags));
nonEscapePosDff = cell(size(folders,2),length(timeLags));
nonEscapeNegDff = cell(size(folders,2),length(timeLags));

%for each assay, get training + testing data for behavior
for mouseNum = 1:length(folders)

    clear escapeFrameMS
    
    dff_neg = DataAll{1}{mouseNum,assayNum}(escapeCellNegId{mouseNum},:);
    dff_pos = DataAll{1}{mouseNum,assayNum}(escapeCellPosId{mouseNum},:);
    
    %z-score cell activity across the session length
    for cellNum = 1:size(dff_neg,1)
        dff_neg(cellNum,:) = zscore(dff_neg(cellNum,:));
    end
    for cellNum = 1:size(dff_pos,1)
        dff_pos(cellNum,:) = zscore(dff_pos(cellNum,:));
    end
    
    if usePCA==1
       X = bsxfun(@minus,dff_neg',mean(dff_neg'));
       [coeff,score,latent,~,explained] = pca(X);
       temp = cumsum(explained); temp = min(find(temp > minVarAcc));
       dff_neg = score(:,1:temp)';
       
       X = bsxfun(@minus,dff_pos',mean(dff_pos'));
       [coeff,score,latent,~,explained] = pca(X);
       temp = cumsum(explained); temp = min(find(temp > minVarAcc));       
       dff_pos = score(:,1:temp)';       
    end
    
    
    %remove behavior if max lag is negative time
    escapeIndices = DataAll{3}{mouseNum,assayNum}(behavNumber,:); %get the escape variable.
    escapeIndices(1) = 0; escapeIndices(end) = 0;
    temp = diff(escapeIndices); escapeFrameMS(:,1) = find(temp==1); escapeFrameMS(:,2) = find(temp==-1);
    clearvars temp
    
    escapeFrameMS(find(escapeFrameMS(:,1) < (max(timeLags)+sampleDur+1)),:) = []; %sampleDur is the amount of data for predicting each escape
    escapeFrameMS(find(escapeFrameMS(:,2) > length(escapeIndices)-sampleDur),:) = []; 
    
    %remove behaviors that happen within 10 seconds of other behaviors
        for behavNum = 1:size(escapeFrameMS,1)-1
                escapeDiffRemoveIdx = diff(escapeFrameMS(:,1));
                escapeDiffRemoveIdx = find(escapeDiffRemoveIdx < (7.5*10)); escapeDiffRemoveIdx = escapeDiffRemoveIdx + 1;
                escapeFrameMS(escapeDiffRemoveIdx,:) = [];
        end
        
        numBehavs(mouseNum) = size(escapeFrameMS,1);

    %now do for each time lag from behavior onset
    for lagNum = 1:length(timeLags)
        %compile all neural data for a given time lag
        for escapeNum = 1:size(escapeFrameMS,1)
            escapePosDff{mouseNum,lagNum} = [escapePosDff{mouseNum,lagNum}, dff_pos(:,[escapeFrameMS(escapeNum,1)-timeLags(lagNum):(escapeFrameMS(escapeNum,1)-timeLags(lagNum)+sampleDur)])];
            escapeNegDff{mouseNum,lagNum} = [escapeNegDff{mouseNum,lagNum}, dff_neg(:,[escapeFrameMS(escapeNum,1)-timeLags(lagNum):(escapeFrameMS(escapeNum,1)-timeLags(lagNum)+sampleDur)])];
        end
    end 
end

%now we need a same-dimensioned cell array with non-escape data.
for mouseNum = 1:length(folders)

    clear escapeFrameMS
    
    dff_neg = DataAll{1}{mouseNum,assayNum}(escapeCellNegId{mouseNum},:);
    dff_pos = DataAll{1}{mouseNum,assayNum}(escapeCellPosId{mouseNum},:);

    %zscore
    for cellNum = 1:size(dff_neg,1)
        dff_neg(cellNum,:) = zscore(dff_neg(cellNum,:));
    end    
    for cellNum = 1:size(dff_pos,1)
        dff_pos(cellNum,:) = zscore(dff_pos(cellNum,:));
    end
    
    if usePCA==1
       X = bsxfun(@minus,dff_neg',mean(dff_neg'));
       [coeff,score,latent,~,explained] = pca(X);
       temp = cumsum(explained); temp = min(find(temp > minVarAcc));
       dff_neg = score(:,1:temp)';
       
       X = bsxfun(@minus,dff_pos',mean(dff_pos'));
       [coeff,score,latent,~,explained] = pca(X);
       temp = cumsum(explained); temp = min(find(temp > minVarAcc));       
       dff_pos = score(:,1:temp)';       
    end    
    
    escapeIndices = DataAll{3}{mouseNum,assayNum}(behavNumber,:); %get the escape variable.
    escapeIndices(1) = 0; escapeIndices(end) = 0;    
    temp = diff(escapeIndices); escapeFrameMS(:,1) = find(temp==1); escapeFrameMS(:,2) = find(temp==-1);
    clearvars temp
    
    %REMEMBER -- don't want to remove escapes that happen within 10 seconds of others -- bc don't want to label them non-escape later
    %make array of possible non-escape indices to sample ('1' means
    %nonescape or non-preescape)
    nonEscapeIndices = ones(1,length(dff_neg));
    
    escapeFrameMS(:,1) = escapeFrameMS(:,1)-max(timeLags);
    escapeFrameMS(find(escapeFrameMS(:,1) < 0),1) = 1;
    
    for escapeNum = 1:size(escapeFrameMS,1)
        nonEscapeIndices([escapeFrameMS(escapeNum,1):escapeFrameMS(escapeNum,2)]) = 0;
    end
    
    indicesToSelect = find(nonEscapeIndices);
    
    %resample this for each time lag
    for lagNum = 1:length(timeLags)
        randIdx = randperm(length(indicesToSelect),size(escapeNegDff{mouseNum,lagNum},2));
        nonEscapeNegDff{mouseNum,lagNum} = dff_neg(:,indicesToSelect(randIdx)); %PS fixed this        
    end    
    for lagNum = 1:length(timeLags)
        randIdx = randperm(length(indicesToSelect),size(escapePosDff{mouseNum,lagNum},2));
        nonEscapePosDff{mouseNum,lagNum} = dff_pos(:,indicesToSelect(randIdx)); %PS fixed this        
    end        
end

%now we have 'escapePosDff' and 'nonEscapePosDff' and vice versa.  Must segment into 5 parts and
%do 5-fold training/testing

for mouseNum = 1:length(folders)
    
    if numBehavs(mouseNum) < 5 %make sure there are at least x escape attempts per mouse
        percCorrectPosMean(mouseNum,1:length(timeLags)) = nan;
        continue
    end

    for lagNum = 1:length(timeLags)
        ['Working on mouse number ', num2str(mouseNum), ', lag number ', num2str(lagNum)]
        clearvars foldIdx
        %find indices for each fold
        lengthTemp = size(escapePosDff{mouseNum,lagNum},2);
        increment = floor(lengthTemp ./ 5);
        foldIdx(1,:) = [1:increment];
        foldIdx(2,:) = [increment.*1+1:increment.*2];
        foldIdx(3,:) = [increment.*2+1:increment.*3];
        foldIdx(4,:) = [increment.*3+1:increment.*4];
        foldIdx(5,:) = [increment.*4+1:increment.*5];
       
        %for each foldNum, remove that group of indices from training and
        %test on it.
            for foldNum = 1:5
                trainingEscData = escapePosDff;
                if foldNum==1
                    trainingEscData{mouseNum,lagNum}(:,foldIdx(foldNum,1):foldIdx(foldNum,end)+26) = []; %remove specific fold indices
                end
                if foldNum==2 | foldNum==3 | foldNum==4
                    trainingEscData{mouseNum,lagNum}(:,foldIdx(foldNum,1)-13:foldIdx(foldNum,end)+13) = []; %remove specific fold indices
                end
                if foldNum==5
                    trainingEscData{mouseNum,lagNum}(:,foldIdx(foldNum,1)-26:foldIdx(foldNum,end)) = []; %remove specific fold indices
                end
               
                trainingNonEscData = nonEscapePosDff;
                if foldNum==1
                    trainingNonEscData{mouseNum,lagNum}(:,foldIdx(foldNum,1):foldIdx(foldNum,end)+26) = []; %remove specific fold indices
                end
                if foldNum==2 | foldNum==3 | foldNum==4
                    trainingNonEscData{mouseNum,lagNum}(:,foldIdx(foldNum,1)-13:foldIdx(foldNum,end)+13) = []; %remove specific fold indices
                end
                if foldNum==5
                    trainingNonEscData{mouseNum,lagNum}(:,foldIdx(foldNum,1)-26:foldIdx(foldNum,end)) = []; %remove specific fold indices
                end

                %and get the testing data for escape and non-escape moments
                testingEscData = escapePosDff; 
                testingEscData{mouseNum,lagNum} = testingEscData{mouseNum,lagNum}(:,foldIdx(foldNum,:)); %remove specific fold indices

                testingNonEscData = nonEscapePosDff; 
                testingNonEscData{mouseNum,lagNum} = testingNonEscData{mouseNum,lagNum}(:,foldIdx(foldNum,:)); %remove specific fold indices

                
                trainingDffAll = [trainingNonEscData{mouseNum,lagNum}, trainingEscData{mouseNum,lagNum}];
                trainingLabelsAll = [zeros(1,size(trainingEscData{mouseNum,lagNum},2))+1, ones(1,size(trainingEscData{mouseNum,lagNum},2))+1];

                testingDffAll = [testingNonEscData{mouseNum,lagNum},testingEscData{mouseNum,lagNum}];
                testingLabelsAll = [zeros(1,size(testingEscData{mouseNum,lagNum},2))+1,ones(1,size(testingEscData{mouseNum,lagNum},2))+1];
                
                B = mnrfit(trainingDffAll',categorical(trainingLabelsAll'));
                pihat = mnrval(B,testingDffAll');
                
                for i = 1:size(pihat,1)
                    [val,loc(i)] = max(pihat(i,:));
                end
                predictedBehav = loc;
                
                %and compare with testing labels to generate a % same
                for i = 1:length(testingLabelsAll)
                    if testingLabelsAll(i) == predictedBehav(i)
                       correct(i) = 1;
                    else
                        correct(i) = 0;
                    end
                end
                percentCorrect(foldNum) = sum(correct) ./ length(correct);
                clearvars correct
            end
          percCorrectPosMean(mouseNum,lagNum) = mean(percentCorrect);
          clearvars percentCorrect loc predictedBehav sum correct
    end     
end

%% Plot the results

        yBound = [.4 1];
        idxPlot = [1,3,5,7,9];
        accPosMean = nanmean(percCorrectPosMean);
        accPosSE = nanstd(percCorrectPosMean) ./ sqrt(size(percCorrectPosMean,1));
        
        %check significance above .5 (chance)
        tempPos = percCorrectPosMean - .5; %make so chance is at '0'.
        for lagNum = 1:size(percCorrectPosMean,2)
            [h,p_Pos(lagNum)] = ttest(tempPos(:,lagNum));
        end
        
        B = subplot(1,1,1);
        bar(accPosMean,'LineWidth',1'); hold on;
        ylim(yBound)
        xlim([.5 length(accPosMean)+.5])
        errorbar(accPosMean,accPosSE,'LineStyle','none','Color','k'); box off;
        title('5-fold logistic regression - predict behavior shock ext.')
        
        for lagNum = 1:size(percCorrectPosMean,2)
            if p_Pos(lagNum) < .05 && p_Pos(lagNum) > .01
                text(lagNum,.9, ['*',' p=',num2str(round(p_Pos(lagNum),4))])
            elseif p_Pos(lagNum) < .01 && p_Pos(lagNum) > .001
                text(lagNum,.9, ['**',' p=',num2str(round(p_Pos(lagNum),4))])
            elseif p_Pos(lagNum) < .001
                text(lagNum,.9, ['***',' p=',num2str(round(p_Pos(lagNum),4))])
            end
        end
        A = plot([0 length(accPosMean)+.5],[.5 .5],':','Color','r'); A.LineWidth = 2; B.LineWidth = 1;
        ylabel('frac. accurate')
        B.XTick = [1:length(accPosMean)];
        B.XTickLabel = round(timeLags ./ 7.5);
        box off;
        xlabel('time lag of -2sec training data from behav onset')


