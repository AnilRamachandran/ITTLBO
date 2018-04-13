function testMode = read_classfi_dataasets(testMode)



[parentdir,~,~] = fileparts(pwd);
dataDir =(strcat(parentdir,'\text data'));
dirData = dir(dataDir);      
dirIndex = [dirData.isdir];  
subDirs = {dirData(dirIndex).name};  
validIndex = ~ismember(subDirs,{'.','..'}); 
for iDir = find(validIndex)                  
    nextDir = fullfile(dataDir,subDirs{iDir}); 
    fileList = dir(fullfile(nextDir,'*.txt'));
    fileNames = {fileList.name};
    train_x = [];train_y = [];test_x = [];test_y = [];
    for ii =1:length(fileNames)
        file= fullfile(nextDir,fileNames{ii});
        data = load(file);
        ri = randperm(size(data,1));
        idxtrain = ri(1:floor(0.7*length(ri)));
        idxtest = ri(floor(0.7*length(ri))+1:end);
        train_x = [train_x;data(idxtrain,1:end-1)];
        train_y = [train_y;data(idxtrain,end)];
        test_x = [test_x;data(idxtest,1:end-1)];
        test_y = [test_y;data(idxtest,end)];
    end 
    testMode.xxtrain{iDir-2} = train_x;
    testMode.xxtest{iDir-2} = test_x;
    testMode.yytrain{iDir-2} = train_y;
    testMode.yytest{iDir-2} = test_y;
end





