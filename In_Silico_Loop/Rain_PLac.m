%% This code should put under AMIGO2 EXAMPLE. I.e. ..\AMIGO2R2016d\Examples\In_Silico_Loop
%% Rain181115 Initial setting
clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899]; % The parameter vector is set
% to the best estimates for MIP,r
%% Rain181115 Off-line OED, sampling time: 2.5, 5, 7.5, 10
Sampling_time_off_line = [2.5, 5, 7.5, 10];
tic
% Sampling_time_off_line = [2.5 5];
errOff = [];
pOff = [];
numLoops = 3;
numExperiments = 1;
resultBaseOff = cell(length(Sampling_time_off_line),1);
parfor i = 1:length(Sampling_time_off_line)
    resultBaseOff{i} = strcat('RainOffSt',strrep(num2str(...
        Sampling_time_off_line(i)),'.','_'));
    run_in_silico_experiment_parfor_Optimised(resultBaseOff{i},numLoops,...
        numExperiments,Sampling_time_off_line(i));
end
for i = 1:length(Sampling_time_off_line)
    load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
    int2str(numExperiments)));
    pOff(i,:) = best_global_theta;
end
errOff = abs(log2(pOff./true_par));
errOffavg = sum(errOff,2)/length(true_par);
bar(Sampling_time_off_line,Sampling_time_off_line)
toc
%% Rain181115 On-line OED, sampling time: 2.5, 5, 7.5, 10
% Note that the code is only for numEcperiments = 3 and numLoops = 3. This
% part should be modified
tic
Sampling_time_on_line = [2.5, 5, 7.5, 10];
errON = [];
pOn = [];
numLoops = 3;
numExperiments = 3;
resultBaseOn = cell(length(Sampling_time_off_line),1);
% pOnLoop1 = ones(length(Sampling_time_on_line),length(true_par));
% pOnLoop2 = ones(length(Sampling_time_on_line),length(true_par));
% pOnLoop3 = ones(length(Sampling_time_on_line),length(true_par));
pOnLoop1 = [];
pOnLoop2 = [];
pOnLoop3 = [];
parfor i = 1:length(Sampling_time_on_line)
    resultBaseOn{i} = strcat('RainOnSt',strrep(num2str(...
        Sampling_time_on_line(i)),'.','_'));
    run_in_silico_experiment_parfor_Optimised(resultBaseOn{i},numLoops,...
        numExperiments,Sampling_time_on_line(i));
end
for i = 1:length(Sampling_time_on_line)
    load(strcat(resultBaseOn{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments-2)));
    pOnLoop1(i,:) = best_global_theta;
    load(strcat(resultBaseOn{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments-1)));
    pOnLoop2(i,:) = best_global_theta;
    load(strcat(resultBaseOn{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments)));
    pOnLoop3(i,:) = best_global_theta;
end
errOnLoop1 = abs(log2(pOnLoop1./true_par));
errOnLoop2 = abs(log2(pOnLoop2./true_par));
errOnLoop3 = abs(log2(pOnLoop3./true_par));
errOnavg = (sum(errOnLoop1,2)/length(true_par)+...
    sum(errOnLoop2,2)/length(true_par)+sum(errOnLoop3,2)/length(true_par))/...
    numExperiments;
bar(Sampling_time_on_line,errOnavg)
toc