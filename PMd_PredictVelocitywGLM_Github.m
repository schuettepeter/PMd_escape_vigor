%% THIS IS THE SCRIPT USED TO PREDICT VELOCITY, USING A GLM, AS DESCRIBED IN THE MANUSCRIPT
%"Dorsal premammillary hypothalamic projection to periaqueductal gray
%controls escape vigor from innate and conditioned threats"

clear all

%ENTER FILEPATHS TO DATA FROM EACH RECORDING SESSION:
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
    
usePCA = 1; %choose to operate on dimensionality-reduced principal components.
minVarAcc = 80; %use the top PCs that account for this amount of variance.

for mouseNum = 1:size(folders,1)
    for assayNum = 1:size(folders,2)
        
        cd(folders{mouseNum,assayNum})
        load('output_CNMF-E.mat','neuron');
        load('BehaviorMS.mat','escapeIndicesMS')
        load('Tracking.mat')
        load('good_neurons.mat')
        %mousePos = Tracking.mouseVelMS(:,1);
        mouseActualPos = Tracking.mouse_positionMS(:,1);
        mousePos = diff(mouseActualPos); %try using velocity, rather than speed
        
        sig = neuron.C_raw(find(good_neurons),:);
        
        if length(mousePos) < length(sig)
            sig = sig(:,1:end-1);
        end
        if length(mousePos) > length(sig)
            mousePos = mousePos(1:end-1);
            mouseActualPos = mouseActualPos(1:end-1);            
        end
        while length(escapeIndicesMS) < length(sig)
            escapeIndicesMS = escapeIndicesMS(:,escapeIndicesMS(end));
        end
        if length(mousePos) > length(sig)
            mousePos = mousePos(1:end-1);
            mouseActualPos = mouseActualPos(1:end-1);
        end
        
        
        if usePCA==1
           X = bsxfun(@minus,sig',mean(sig'));
           [coeff,score,latent,~,explained] = pca(X);
           temp = cumsum(explained); temp = min(find(temp > minVarAcc));
           sig = score(:,1:temp)';
        end
        
        mousePos = smoothdata(mousePos); %smooth the velocity

    if 0==1    
    %truncate so only periods where mouse is closer to threat
    idxToDel = find(mouseActualPos > 425);
    sig(:,idxToDel) = [];
    mousePos(idxToDel) = [];
    escapeIndicesMS(idxToDel) = [];
    end
    
    %split data 
    dataSegLength = floor(60 .* 7.5);  
    inbetweenSegs = floor(10 .* 7.5);
    iterTotal = floor(length(sig) ./ dataSegLength);
    %iterTotal = floor(iterTotal ./ 2);        
    
    trainIdx = zeros(1,length(sig)); 
    testIdx = zeros(1,length(sig));   
        for iterNum = 1:iterTotal       
           if bitget(iterNum,1) %odd
                trainIdx(((iterNum-1).*dataSegLength)+1:(iterNum .* dataSegLength-inbetweenSegs)) = 1;
                testIdx(((iterNum-1).*dataSegLength)+1:(iterNum .* dataSegLength-inbetweenSegs)) = 0;            
           else %even
                trainIdx(((iterNum-1).*dataSegLength)+1:(iterNum .* dataSegLength-inbetweenSegs)) = 0;
                testIdx(((iterNum-1).*dataSegLength)+1:(iterNum .* dataSegLength-inbetweenSegs)) = 1;            
           end 
        end
        
     sigTrain = sig(:,find(trainIdx));
     posTrain = mousePos(find(trainIdx));
     sigTest = sig(:,find(testIdx));
     posTest = mousePos(find(testIdx));
        
     mdl = fitglm(sigTrain',posTrain');
     posPredict = predict(mdl,sigTest');
     
     error = posPredict-posTest;
     
     MSE(mouseNum,assayNum) = nansum(error.^2) ./ length(sigTest);
     
     posPredTestAll{mouseNum,assayNum} = [posPredict,posTest];

    end
end

%% plot the MSE for prediction, control vs. threat
%first convert to cm
[p_rat,h] = signrank(MSE(:,1),MSE(:,2));
[p_shk,h] = signrank(MSE(:,3),MSE(:,4));

figure(135)
bar(nanmean(MSE)); hold on;
errorbar(nanmean(MSE),(nanstd(MSE))./sqrt(9), 'LineStyle','none','Color','k')
box off;
title('GLM MSE mouse velocity: toy=1,rat=2,pre-shk=3,post-shk=4')
text(1.5,7,num2str(p_rat));
text(3.5,7,num2str(p_shk));
ylim([0 8])