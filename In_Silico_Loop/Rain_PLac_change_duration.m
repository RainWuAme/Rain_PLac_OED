%% This code should put under AMIGO2 EXAMPLE. I.e. ..\AMIGO2R2016d\Examples\In_Silico_Loop
%% Rain181120 Off-line OED, sampling time: 2.5, 5, 10
clear all
dbstop if error
profile on
Sampling_time_off_line = [5,6];
% 1000/sampling_time and 300/sampling_time -> the result needs to be
% integer.
numLoops = 1; % Number of Sub-experiment. 1 means there is only one 
% experiment. That means off-line OED.
numExperiments = 1; % How many times of experiment
resultBaseOff = cell(length(Sampling_time_off_line),1); % Preallocation

% File name input setting. The input will be, for example, RainOffSt2_5.
for i = 1:length(Sampling_time_off_line)
    resultBaseOff{i} = strcat('RainOffSt',strrep(num2str(...
        Sampling_time_off_line(i)),'.','_'));
end

% OED
parfor i = 1:length(Sampling_time_off_line)
    run_in_silico_experiment_parfor_Optimised_light(resultBaseOff{i},numLoops,...
        numExperiments,Sampling_time_off_line(i));
end
p = profile('info');
%% Rain181115 On-line OED, sampling time: 2.5, 5,  10
%{
clear all, close all, clc
dbstop if error
profile on
% Sampling_time_on_line = [2.5, 5, 10];
Sampling_time_on_line = [5 10];
numLoops = 3;
% numExperiments = 3;
numExperiments = 1;
resultBaseOn = cell(length(Sampling_time_on_line),1);

% File name input setting. The input will be, for example, RainOnSt2_5.
for i = 1:length(Sampling_time_on_line)
    resultBaseOn{i} = strcat('RainOnSt',strrep(num2str(...
        Sampling_time_on_line(i)),'.','_'));
end

% OED
parfor i = 1:length(Sampling_time_on_line)
    run_in_silico_experiment_parfor_Optimised(resultBaseOn{i},numLoops,...
        numExperiments,Sampling_time_on_line(i));
end
profsave
%% Plot result off-line (Edinburgh)
clear all, clc, close all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
cd('D:\AMIGO2R2016d\Examples\In_Silico_Loop\Edinburgh_lab_result')
Sampling_time_off_line = [2.5, 5, 10];
numLoops = 1; % Number of Sub-experiment. 1 means there is only one 
% experiment. That means off-line OED.
numExperiments = 3; % How many times of experiment
pOff = [];
resultBaseOff = cell(length(Sampling_time_off_line),1);
for i = 1:length(Sampling_time_off_line)
    resultBaseOff{i} = strcat('RainOffSt',strrep(num2str(...
    Sampling_time_off_line(i)),'.','_'));
end
for i = 1:length(Sampling_time_off_line)
    for j = 1:numExperiments
        load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(j)))
        pOff(i,:,j) = best_global_theta';
    end
end
% pOff contains "numExperiments" 3D matrix. Each of it represent the result of
% experimtent 1, experiment 2 and so on. Each of the cell contains a
% matrix which the raws represent different sampling time and the column
% represent each of the parameters.

errOff = abs(log2(pOff./true_par));
errOffavg = sum(sum(errOff,2)/length(true_par),3)/numExperiments;

figure(1)
[colormap] = cbrewer('qual','Set1',3); % Color brewer
hb = bar(1:length(Sampling_time_off_line),errOffavg,'FaceColor',...
    colormap(2,:)); hold on
% errorbar(hb.XData,hb.YData,[1.96.*std(reshape(errOff,...
%     length(Sampling_time_off_line),[]),0,2)],'k.')
set(gca,'xticklabel',Sampling_time_off_line)
xlabel('Sampling time','FontSize',10,'Interpreter','Latex')
ylabel('$\bar{\epsilon}$','FontSize',10,'Interpreter','Latex')
title('Off-Line OED, 3 exp, 1 subexp, bar plot (Edinburgh)')
% ylim([-1.5,2])
ylim([0,0.5])
% legend('Average relative error','95% interval')
legend('Average relative error')
set(gcf, 'PaperPosition', [0 0 13 13]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [13 13])
saveas(gca,'Edinburgh_Off_line_OED_3_sampling_time.pdf');


figure(2)
boxplot(permute(errOff(1,:,:),[3,2,1])); hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
% title('Off-Line OED, St2.5, 3 exp, 1 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

figure(3)
boxplot(permute(errOff(2,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
% title('Off-Line OED, St5, 3 exp, 1 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

figure(4)
boxplot(permute(errOff(3,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
% title('Off-Line OED, St10, 3 exp, 1 subexp, parameters box plot (Edinburgh)')
ylim([0,3])
%% Plot result on-line (Edinburgh)
clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
cd('D:\AMIGO2R2016d\Examples\In_Silico_Loop\Edinburgh_lab_result')
Sampling_time_on_line = [2.5, 5, 10];
numLoops = 3; % Subexperiment
numExperiments = 3; % Do 3 times experiment
resultBaseOn = cell(length(Sampling_time_on_line),1);
pOn = [];
for i = 1:length(Sampling_time_on_line)
    resultBaseOn{i} = strcat('RainOnSt',strrep(num2str(...
    Sampling_time_on_line(i)),'.','_'));
end
for i = 1:length(Sampling_time_on_line)
    for j = 1:numExperiments
        load(strcat(resultBaseOn{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(j)))
        pOn(i,:,j) = best_global_theta';
    end
end
% The row of errOffLoopo1~3 is the sampling time. The column is the best
% parameters.
errOn = abs(log2(pOn./true_par));
errOnavg = sum(sum(errOn,2)/length(true_par),3)/numExperiments;

figure(1)
[colormap] = cbrewer('qual','Set1',3);
hb = bar(1:length(Sampling_time_on_line),errOnavg,'FaceColor',colormap(2,:));
hold on
% errorbar(hb.XData,hb.YData,[1.96.*std(reshape(errOn,...
%     length(Sampling_time_on_line),[]),0,2)],'k.')
set(gca,'xticklabel',Sampling_time_on_line)
xlabel('Sampling time','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','Interpreter','Latex')
% title('On-Line OED, 3 exp, 3subexp, error bar plot (Edinburgh)')
ylim([0,0.5])
legend('Average relative error')
set(gcf, 'PaperPosition', [0 0 13 13]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [13 13])
saveas(gca,'Edinburgh_On_line_OED_3_sampling_time.pdf');

figure(2)
boxplot(permute(errOn(1,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St2.5, 3 exp, 3 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

figure(3)
boxplot(permute(errOn(2,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St5, 3 exp, 3 subexp, parameters box plot (Edinburgh)')
ylim([0,3])

figure(4)
boxplot(permute(errOn(3,:,:),[3,2,1])), hold on
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
cd('D:\AMIGO2R2016d\Examples\In_Silico_Loop\20181125')
Sampling_time_off_line = [2.5, 5, 10];
numLoops = 1; % Number of Sub-experiment. 1 means there is only one 
% experiment. That means off-line OED.
numExperiments = 3; % How many times of experiment
pOff = [];
resultBaseOff = cell(length(Sampling_time_off_line),1);
for i = 1:length(Sampling_time_off_line)
    resultBaseOff{i} = strcat('RainOffSt',strrep(num2str(...
    Sampling_time_off_line(i)),'.','_'));
end
for i = 1:length(Sampling_time_off_line)
    for j = 1:numExperiments
        load(strcat(resultBaseOff{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(j)))
        pOff(i,:,j) = best_global_theta';
    end
end
% pOff contains "numExperiments" 3D matrix. Each of it represent the result of
% experimtent 1, experiment 2 and so on. Each of the cell contains a
% matrix which the raws represent different sampling time and the column
% represent each of the parameters.

errOff = abs(log2(pOff./true_par));
errOffavg = sum(sum(errOff,2)/length(true_par),3)/numExperiments;

figure(1)
[colormap] = cbrewer('qual','Set1',3);
hb = bar(1:length(Sampling_time_off_line),errOffavg,'FaceColor',...
    colormap(2,:)); hold on
% errorbar(hb.XData,hb.YData,[1.96.*std(reshape(errOff,...
%     length(Sampling_time_off_line),[]),0,2)],'k.')
set(gca,'xticklabel',Sampling_time_off_line)
xlabel('Sampling time','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','Interpreter','Latex')
title('Off-Line OED, 3 exp, 1 subexp, error bar plot (Rain)')
ylim([0,0.5])
legend('Average relative error')

figure(2)
boxplot(permute(errOff(1,:,:),[3,2,1])); hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St2.5, 3 exp, 1 subexp, parameters box plot (Rain)')
ylim([0,3])

figure(3)
boxplot(permute(errOff(2,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St5, 3 exp, 1 subexp, parameters box plot (Rain)')
ylim([0,3])

figure(4)
boxplot(permute(errOff(3,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('Off-Line OED, St10, 3 exp, 1 subexp, parameters box plot (Rain)')
ylim([0,3])
%% Plot result on-line (Rain)
clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
cd('D:\AMIGO2R2016d\Examples\In_Silico_Loop\20181125')
Sampling_time_on_line = [2.5, 5, 10];
numLoops = 3; % Subexperiment
numExperiments = 3; % Do 3 times experiment
resultBaseOn = cell(length(Sampling_time_on_line),1);
pOn = [];
for i = 1:length(Sampling_time_on_line)
    resultBaseOn{i} = strcat('RainOnSt',strrep(num2str(...
    Sampling_time_on_line(i)),'.','_'));
end
for i = 1:length(Sampling_time_on_line)
    for j = 1:numExperiments
        load(strcat(resultBaseOn{i},'-OptstepseSS-',int2str(numLoops),'_loops-',...
        int2str(j)))
        pOn(i,:,j) = best_global_theta';
    end
end
% The row of errOffLoopo1~3 is the sampling time. The column is the best
% parameters.
errOn = abs(log2(pOn./true_par));
errOnavg = sum(sum(errOn,2)/length(true_par),3)/numExperiments;

figure(1)
[colormap] = cbrewer('qual','Set1',3);
hb = bar(1:length(Sampling_time_on_line),errOnavg,'FaceColor',colormap(2,:));
hold on
% errorbar(hb.XData,hb.YData,[1.96.*std(reshape(errOn,...
%     length(Sampling_time_on_line),[]),0,2)],'k.')
set(gca,'xticklabel',Sampling_time_on_line)
xlabel('Sampling time','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','Interpreter','Latex')
% title('On-Line OED, 3 exp, 3subexp, error bar plot (Edinburgh)')
ylim([0,0.5])
legend('Average relative error')
set(gcf, 'PaperPosition', [0 0 13 13]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [13 13])
saveas(gca,'Edinburgh_On_line_OED_3_sampling_time.pdf');

figure(2)
boxplot(permute(errOn(1,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','Interpreter','Latex')
title('On-Line OED, St2.5, 3 exp, 3 subexp, parameters box plot (Rain)')
ylim([0,3])

figure(3)
boxplot(permute(errOn(2,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St5, 3 exp, 3 subexp, parameters box plot (Rain)')
ylim([0,3])

figure(4)
boxplot(permute(errOn(3,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St10, 3 exp, 3 subexp, parameters box plot (Rain)')
ylim([0,3])
%% Plot OID result - examples
close all
% Sampling time 2.5, experiment 1, on-line 
[colormap] = cbrewer('qual','Set1',3); % Color brewer
figure(5)
load('RainOffSt2_5-OptstepseSS-1_loops-1.mat')
IPTG = exps.u{1,1};
t_switch = exps.t_con{1,1};
t_plot = repelem(t_switch,2)/60;
t_plot(1) = [];
t_plot(end) = [];
t = oed_results{1,1}.sim.tsim{1}/60;
yyaxis('left')
plot(t,inputs.exps.exp_data{1,1},'.','Color',colormap(3,:)), hold on
plot(t,pe_results{1,1}.sim.obs{1,1},'-','Color',colormap(3,:),'LineWidth',2)
ax = gca;
ax.YColor = colormap(3,:);
ylim([-100 1200])
xlabel('Time (hours)')
ylabel('Citrine ($10^3$ molecules)','Interpreter','Latex','Color',colormap(3,:))
yyaxis('right')
IPTG_plot = repelem(IPTG,2);
plot(t_plot,IPTG_plot,'Color',colormap(1,:),'LineWidth',2)
ylim([-100 1200])
ylabel('IPTG ($10^3\mu M$ )','Color',colormap(1,:),'Interpreter','Latex')
title('OID result - sampling time 2.5, experiment 1, off-line')
% Note that the model structure is M_{3D} in the paper.

% Sampling time 5, experiment 2, on-line
figure(6)
load('RainOnSt5-OptstepseSS-3_loops-3')
IPTG = exps.u{1,1};
t_switch = exps.t_con{1,1};
t_plot = repelem(t_switch,2)/60;
t_plot(1) = [];
t_plot(end) = [];
t = inputs.exps.t_s{1,1}/60;
yyaxis('left')
plot(t,inputs.exps.exp_data{1,1},'.','Color',colormap(3,:)), hold on
plot(t,pe_results{1,3}.sim.obs{1,1},'-','Color',colormap(3,:),'LineWidth',2)
ax = gca;
ax.YColor = colormap(3,:);
ylim([-100 1200])
xlabel('Time (hours)')
ylabel('Citrine ($10^3$ molecules)','Interpreter','Latex','Color',colormap(3,:))
yyaxis('right')
IPTG_plot = repelem(IPTG,2);
plot(t_plot,IPTG_plot,'Color',colormap(1,:),'LineWidth',2)
ylim([-100 1200])
ylabel('IPTG ($10^3\mu M$ )','Color',colormap(1,:),'Interpreter','Latex')
title('OID result - sampling time 5, experiment 3, on-line')
%}
%% Test spmd
clear all
% dbstop if error
profile on
Sampling_time_off_line = 5;
% 1000/sampling_time and 300/sampling_time -> the result needs to be
% integer.
numLoops = 1; % Number of Sub-experiment. 1 means there is only one 
% experiment. That means off-line OED.
numExperiments = 1; % How many times of experiment
resultBaseOff = cell(length(Sampling_time_off_line),1); % Preallocation

% File name input setting. The input will be, for example, RainOffSt2_5.
for i = 1:length(Sampling_time_off_line)

end

% OED
try
    parpool(2)
catch
    disp('Parpool already on')
end

spmd
    resultBaseOff = strcat('RainOffSt',strrep(num2str(Sampling_time_off_line),'.','_'),'Index',num2str(labindex));
    run_in_silico_experiment_parfor_Optimised_light(char(resultBaseOff(labindex)),numLoops,numExperiments,Sampling_time_off_line);
end

p = profile('info');