clear all;
close all;
mixingCoff = [];
str = 'elastic_net';
if isempty(str)
    str = 'synthetic';
end
for outerLoop = 1 : 10
    BO_TL_PES_loop;
    noTransferY(:,outerLoop) = (dataTarget.y*sqrt(dataTarget.var_y))+dataTarget.mean_y;
    T2TransferY(:,outerLoop) = (dataTargetT2.y*sqrt(dataTargetT2.var_y))+dataTargetT2.mean_y;
    EnvTransferY(:,outerLoop) = (dataTarget_MergedSrc.y*sqrt(dataTarget_MergedSrc.var_y))+dataTarget_MergedSrc.mean_y;
    YogTransferY(:,outerLoop) = (dataTargetYog.y0*sqrt(dataTargetYog.var_y))+dataTargetYog.mean_y;
    BrdTransferY(:,outerLoop) = (dataTargetBrd.y0*sqrt(dataTargetBrd.var_y))+dataTargetBrd.mean_y;
    
    mixingCoff(:,:,outerLoop) = [percentage_matrix];
    
    clearvars -except outerLoop noTransferY T2TransferY EnvTransferY YogTransferY BrdTransferY str mixingCoff 
end
%% Performace Graph
itersRun = size(T2TransferY,2);
h=figure;
x=1:size(T2TransferY,1);
y = randn(4,length(T2TransferY),itersRun);
y(1,:,:)=cummin(T2TransferY);y(2,:,:)=cummin(EnvTransferY);y(3,:,:)=cummin(YogTransferY);y(4,:,:)=cummin(BrdTransferY);
y(5,:,:)=cummin(noTransferY);
lineProps.col{1} = 'b';
lineProps.col{2} = 'm';
lineProps.col{3} = 'g';
lineProps.col{4} = 'r';
lineProps.col{5} = 'k';
H = mseb(x,mean(y,3), std(y,[],3)/sqrt(outerLoop), lineProps);
legend('Proposed Method','Env-GP','SMBO','SCoT','Generic-BO','Location','NorthEast');
xlabel('Number of  Iterations');ylabel('Minimum Value');
%% Percentage bar graph
q=figure;
bar(mean(mixingCoff,3)','stacked');
legend('Source 1','Source 2','Source 3','Source 4','Target');