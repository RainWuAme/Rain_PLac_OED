function [out] = Rain_PLac_wrapper(StepNum,NumSubExp,ExpIndex)
% Rain_PLac_wrapper(StepNum,NumSubExp,ExpIndex)
% StepNum: Number of step input. This should be a postive integer. Note
% that the duration of the all experiment is 3000 min.
% 
% NumSubExp: An array define the number of sub-experiment (on-line OED). 
% The element in it needs to be the factor of StepNum. E.g. If StepNum = 15 
% then NumSubExp can be the combinaiton of 1, 3, 5, and 15 such as 
% [1,3,5,15] or [3,5] or [1,3,15] and so on.
% 
% ExpIndex: The index of the experiment. 
% 
% This is a monte carlo simulation. If StepNum = 40, NumSubExp = [4], and 
% ExpIndex = [1,2,5] then the result file will be named as
% Rain_StepNum-40_NumSubExp-4_ExpIndex-1,
% Rain_StepNum-40_NumSubExp-4_ExpIndex-2,
% Rain_StepNum-40_NumSubExp-4_ExpIndex-5

true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899]; % The parameter vector is set
% to the best estimates for MIP,r
Sampling_time_on_line = 5;
numExperiments = 1; % Do 1 times experiment
resultBaseOn = cell(length(ExpIndex),length(NumSubExp));
for i = 1:length(ExpIndex)
    for j = 1:length(NumSubExp)
    resultBaseOn{i,j} = strcat('Rain_Step',num2str(StepNum),...
        '_SubExp',num2str(NumSubExp(j)),'_Exp',...
        num2str(ExpIndex(i)));
    Rain_run_in_silico_experiment_parfor_Optimised(resultBaseOn{i,j},...
        NumSubExp(j),numExperiments,Sampling_time_on_line,StepNum);
    end
end
%     run_in_silico_experiment_parfor_Optimised(resultBaseOff{i},numLoops,...
%         numExperiments,Sampling_time_off_line(i));
out = 1;
