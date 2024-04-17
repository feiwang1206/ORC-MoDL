function prepare_input_label(pname_net,fname,n_epoch,lr)
%%
    addpath(genpath('/home/ubuntu/Documents/MATLAB/mfiles'));

    epoch(1)=n_epoch;
    MiniBatchSize(1)=2;
    
    for ii=1:1
        trainfolder{1}=[fname{ii}{1} '/inputsTra'];
        trainfolder{2}=[fname{ii}{2} '/labelsTra'];
        for ii2=1:length(trainfolder)
            volReader = @(x) matRead2(x);
            volLoc = trainfolder{ii2};
            inputsTemp = imageDatastore(volLoc, 'FileExtensions','.mat','ReadFcn',volReader);
            inputs{ii2}=inputsTemp;
        end

        trainfolder{1}=[fname{ii}{1} '/inputsVal'];
        trainfolder{2}=[fname{ii}{2} '/labelsVal'];
        for ii2=1:length(trainfolder)
            volReader = @(x) matRead2(x);
            volLoc = trainfolder{ii2};
            inputsTemp = imageDatastore(volLoc, 'FileExtensions','.mat','ReadFcn',volReader);
            inputsVal{ii2}=inputsTemp;
        end
        
        inputname=inputs{1}.Files(1);
        imageSize{1}=size(matRead2(inputname{1}));
        if length(imageSize{1})<4
            imageSize{1}(4)=1;
        end

        inputname=inputs{2}.Files(1);
        channel=size(matRead2(inputname{1}));
        if length(channel)<4
            channel(4)=1;
        end
        
        ds = combine(inputs{1},inputs{2});
        dsVal = combine(inputsVal{1},inputsVal{2});
        mkdir(pname_net);

        options = trainingOptions('adam', ...
            'InitialLearnRate',lr(ii), ...
            'LearnRateSchedule','piecewise',...
            'LearnRateDropPeriod',10,...
            'LearnRateDropFactor',0.5,...
            'Shuffle','every-epoch',...
            'MaxEpochs',n_epoch, ...
            'VerboseFrequency',floor(length(inputs{1}.Files)/MiniBatchSize(ii)),...
            'MiniBatchSize',2,...
            'Plots','training-progress',...
            'ValidationFrequency',floor(length(inputs{1}.Files)/MiniBatchSize(ii)),...
            'ValidationData',dsVal);
        if exist([pname_net '/net_block.mat'],'file')
            load([pname_net '/net_block.mat'],'net');
            [net,info] = trainNetwork(ds,layerGraph(net),options);
        else
            lgraph = createDiffusionNetwork3d_train(imageSize{1}(1:3),imageSize{1}(4),channel(4));
            [net,info] = trainNetwork(ds,lgraph,options);
        end
        save([pname_net '/net_block.mat'],'net');
        save([pname_net '/info.mat'],'info');
    end

end
