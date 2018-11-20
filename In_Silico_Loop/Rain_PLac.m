%% This code should put under AMIGO2 EXAMPLE. I.e. ..\AMIGO2R2016d\Examples\In_Silico_Loop
%% Rain181115 Initial setting
clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899]; % The parameter vector is set
% to the best estimates for MIP,r
%% Rain181115 Off-line OED, sampling time: 2.5, 5, 7.5, 10
Sampling_time_off_line = [2.5, 5, 10];
% 1000/sampling_time and 300/sampling_time -> the result needs to be
% integer.
dbstop if error
errOff = [];
pOff = [];
numLoops = 1; % Subexperiment
numExperiments = 3; % Do 3 times experiment
resultBaseOff = cell(length(Sampling_time_off_line),1);
for i = 1:length(Sampling_time_off_line)
    resultBaseOff{i} = strcat('RainOffSt',strrep(num2str(...
    Sampling_time_off_line(i)),'.','_'));
end
parfor i = 1:length(Sampling_time_off_line)
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

figure(1)
[colormap] = cbrewer('qual','Set1',3);
bar(1:length(Sampling_time_off_line),errOffavg,'FaceColor',colormap(2,:))
set(gca,'xticklabel',Sampling_time_off_line)
xlabel('Sampling time','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, 3 loops, error bar plot (Rain)')

figure(2)
boxplot(errOff), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, 3 loops, parameters box plot (Rain)')
%% Rain181115 On-line OED, sampling time: 2.5, 5, 7.5, 10
% Note that the code is only for numEcperiments = 3 and numLoops = 3. This
% part should be modified
tic
Sampling_time_on_line = [2.5, 5, 10];
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
%% Plot result off-line
clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
cd('D:\AMIGO2R2016d\Examples\In_Silico_Loop\Edinburgh_lab_result')
Sampling_time_off_line = [2.5, 5, 10];
errOff = [];
pOff = [];
numLoops = 1; % Subexperiment
numExperiments = 3; % Do 3 times experiment
resultBaseOff = cell(length(Sampling_time_off_line),1);
for i = 1:length(Sampling_time_off_line)
    resultBaseOff{i} = strcat('RainOffSt',strrep(num2str(...
    Sampling_time_off_line(i)),'.','_'));
end
for i = 1:length(Sampling_time_off_line)
    load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments-2)));
    pOffLoop1(i,:) = best_global_theta;
    load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments-1)));
    pOffLoop2(i,:) = best_global_theta;
    load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments)));
    pOffLoop3(i,:) = best_global_theta;
end
% The row of errOffLoopo1~3 is the sampling time. The column is the best
% parameters.
errOffLoop1 = abs(log2(pOffLoop1./true_par));
errOffLoop2 = abs(log2(pOffLoop2./true_par));
errOffLoop3 = abs(log2(pOffLoop3./true_par));
errOffavg = (sum(errOffLoop1,2)/length(true_par)+...
    sum(errOffLoop2,2)/length(true_par)+sum(errOffLoop3,2)/length(true_par))/...
    numExperiments;

figure(1)
[colormap] = cbrewer('qual','Set1',3);
hb = bar(1:length(Sampling_time_off_line),errOffavg,'FaceColor',...
    colormap(2,:)); hold on
errorbar(hb.XData,hb.YData,[1.96*std([errOffLoop1(1,:),errOffLoop2(1,:),...
    errOffLoop3(1,:)]),1.96*std([errOffLoop1(2,:),errOffLoop2(2,:),...
    errOffLoop3(2,:)]),1.96*std([errOffLoop1(3,:),errOffLoop2(3,:),...
    errOffLoop3(3,:)])],'k.')
set(gca,'xticklabel',Sampling_time_off_line)
xlabel('Sampling time','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, 3 exp, 1 subexp, error bar plot (Edinburgh)')
ylim([-1.5,2])
plot(1,[errOffLoop1(1,:),errOffLoop2(1,:),errOffLoop3(1,:)],'k*')
plot(2,[errOffLoop1(2,:),errOffLoop2(2,:),errOffLoop3(2,:)],'k*')
plot(3,[errOffLoop1(3,:),errOffLoop2(3,:),errOffLoop3(3,:)],'k*')
legend('Average relative error','95% interval','Relative error')

figure(2)
boxplot([errOffLoop1(1,:); errOffLoop2(1,:); errOffLoop3(1,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St2.5, 3 exp, 1 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

figure(3)
boxplot([errOffLoop1(2,:); errOffLoop2(2,:); errOffLoop3(2,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St5, 3 exp, 1 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

figure(4)
boxplot([errOffLoop1(3,:); errOffLoop2(3,:); errOffLoop3(3,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St10, 3 exp, 1 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

%% Plot result on-line
clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
cd('D:\AMIGO2R2016d\Examples\In_Silico_Loop\Edinburgh_lab_result')
Sampling_time_on_line = [2.5, 5, 10];
numLoops = 3; % Subexperiment
numExperiments = 3; % Do 3 times experiment
resultBaseOn = cell(length(Sampling_time_on_line),1);
for i = 1:length(Sampling_time_on_line)
    resultBaseOn{i} = strcat('RainOnSt',strrep(num2str(...
    Sampling_time_on_line(i)),'.','_'));
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
% The row of errOffLoopo1~3 is the sampling time. The column is the best
% parameters.
errOnLoop1 = abs(log2(pOnLoop1./true_par));
errOnLoop2 = abs(log2(pOnLoop2./true_par));
errOnLoop3 = abs(log2(pOnLoop3./true_par));
errOnavg = (sum(errOnLoop1,2)/length(true_par)+...
    sum(errOnLoop2,2)/length(true_par)+sum(errOnLoop3,2)/length(true_par))/...
    numExperiments;

figure(1)
[colormap] = cbrewer('qual','Set1',3);
hb = bar(1:length(Sampling_time_on_line),errOnavg,'FaceColor',colormap(2,:));
hold on
errorbar(hb.XData,hb.YData,[1.96*std([errOnLoop1(1,:),errOnLoop2(1,:),...
    errOnLoop3(1,:)]),1.96*std([errOnLoop1(2,:),errOnLoop2(2,:),...
    errOnLoop3(2,:)]),1.96*std([errOnLoop1(3,:),errOnLoop2(3,:),...
    errOnLoop3(3,:)])],'k.')
set(gca,'xticklabel',Sampling_time_on_line)
xlabel('Sampling time','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, 3 exp, 3subexp, error bar plot (Edinburgh)')
ylim([-1.5,2])
plot(1,[errOnLoop1(1,:),errOnLoop2(1,:),errOnLoop3(1,:)],'k*')
plot(2,[errOnLoop1(2,:),errOnLoop2(2,:),errOnLoop3(2,:)],'k*')
plot(3,[errOnLoop1(3,:),errOnLoop2(3,:),errOnLoop3(3,:)],'k*')
legend('Average relative error','95% interval','Relative error')

figure(2)
boxplot([errOnLoop1(1,:); errOnLoop2(1,:); errOnLoop3(1,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St2.5, 3 exp, 3 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

figure(3)
boxplot([errOnLoop1(2,:); errOnLoop2(2,:); errOnLoop3(2,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St5, 3 exp, 3 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

figure(4)
boxplot([errOnLoop1(3,:); errOnLoop2(3,:); errOnLoop3(3,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St10, 3 exp, 3 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

%% Plot result off-line (Rain)
clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
Sampling_time_off_line = [2.5, 5, 10];
cd('D:\AMIGO2R2016d\Examples\In_Silico_Loop')
errOff = [];
pOff = [];
numLoops = 1; % Subexperiment
numExperiments = 3; % Do 3 times experiment
resultBaseOff = cell(length(Sampling_time_off_line),1);
for i = 1:length(Sampling_time_off_line)
    resultBaseOff{i} = strcat('RainOffSt',strrep(num2str(...
    Sampling_time_off_line(i)),'.','_'));
end
for i = 1:length(Sampling_time_off_line)
    load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments-2)));
    pOffLoop1(i,:) = best_global_theta;
    load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments-1)));
    pOffLoop2(i,:) = best_global_theta;
    load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(numExperiments)));
    pOffLoop3(i,:) = best_global_theta;
end
% The row of errOffLoopo1~3 is the sampling time. The column is the best
% parameters.
errOffLoop1 = abs(log2(pOffLoop1./true_par));
errOffLoop2 = abs(log2(pOffLoop2./true_par));
errOffLoop3 = abs(log2(pOffLoop3./true_par));
errOffavg = (sum(errOffLoop1,2)/length(true_par)+...
    sum(errOffLoop2,2)/length(true_par)+sum(errOffLoop3,2)/length(true_par))/...
    numExperiments;

figure(1)
[colormap] = cbrewer('qual','Set1',3);
hb = bar(1:length(Sampling_time_off_line),errOffavg,'FaceColor',colormap(2,:));
hold on
errorbar(hb.XData,hb.YData,[1.96*std([errOffLoop1(1,:),errOffLoop2(1,:),...
    errOffLoop3(1,:)]),1.96*std([errOffLoop1(2,:),errOffLoop2(2,:),...
    errOffLoop3(2,:)]),1.96*std([errOffLoop1(3,:),errOffLoop2(3,:),...
    errOffLoop3(3,:)])],'k.')
set(gca,'xticklabel',Sampling_time_off_line)
xlabel('Sampling time','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, 3 exp, 1 subexp, error bar plot (Rain)')
ylim([-1.5,2])
plot(1,[errOffLoop1(1,:),errOffLoop2(1,:),errOffLoop3(1,:)],'k*')
plot(2,[errOffLoop1(2,:),errOffLoop2(2,:),errOffLoop3(2,:)],'k*')
plot(3,[errOffLoop1(3,:),errOffLoop2(3,:),errOffLoop3(3,:)],'k*')
legend('Average relative error','95% interval','Relative error')

figure(2)
boxplot([errOffLoop1(1,:); errOffLoop2(1,:); errOffLoop3(1,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St2.5, 3 exp, 1 subexp, parameters box plot (Rain)')
ylim([0,3])

figure(3)
boxplot([errOffLoop1(2,:); errOffLoop2(2,:); errOffLoop3(2,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St5, 3 exp, 1 subexp, parameters box plot (Rain)')
ylim([0,3])

figure(4)
boxplot([errOffLoop1(3,:); errOffLoop2(3,:); errOffLoop3(3,:)]), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St10, 3 exp, 1 subexp, parameters box plot (Rain)')
ylim([0,3])
