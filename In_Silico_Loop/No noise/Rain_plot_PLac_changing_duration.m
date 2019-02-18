clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
cd('D:\rain\PLac_result\Data_change_duration')

% Initialization
step_num = 40;
Sampling_time_on_line = 5;
numLoops = [1,4,8,10,20]; % Subexperiment
numExperiments = 30; % Do 3 times experiment
resultBaseOn = {};
par = [];
SubExp = {};

% Load data, SubExp -> row: numLoops, column: index of experiments
% par is 3 dimentional matrix -> row: numLoops, column: best theta, z: 30 times
% experiment
for i = 1:length(numLoops)
    for j = 1:numExperiments
        SubExp{i,j} = load(strcat('Rain_Step', int2str(step_num),...
            '_SubExp', int2str(numLoops(i)), '_Exp', int2str(j)));
        par(i,:,j) = SubExp{i,j}.best_global_theta';
    end
end

% The row of errOffLoopo1~3 is the sampling time. The column is the best
% parameters.
errOn = abs(log2(par./true_par));
errOnavg = sum(sum(errOn,2)/length(true_par),3)/numExperiments;

figure(1)
[colormap] = cbrewer('qual','Set1',3);
hb = bar(numLoops,errOnavg,'FaceColor',colormap(2,:));
hold on
% errorbar(hb.XData,hb.YData,[1.96.*std(reshape(errOn,...
%     length(Sampling_time_on_line),[]),0,2)],'k.')
% set(gca,'xticklabel',numLoops)
xlabel('Number of subexperiment','Interpreter','Latex')
ylabel('$\bar{\epsilon}$','Interpreter','Latex')
title('On-Line OED, 30 times exp, Sampling time 5 min')
legend('Average relative error')
set(gcf, 'PaperPosition', [0 0 13 13]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [13 13])
ylim([0,0.15])
saveas(gca,'Edinburgh_On_line_OED_3_sampling_time.pdf');

figure(2)
boxplot(permute(errOn(1,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','Interpreter','Latex')
title('Off-Line OED, St5, 30 exp, 1 subexp, parameters box plot (Rain)')
ylim([0,1])

figure(3)
boxplot(permute(errOn(2,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St5, 30 exp, 4 subexp, parameters box plot (Rain)')
ylim([0,1])

figure(4)
boxplot(permute(errOn(3,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St5, 30 exp, 8 subexp, parameters box plot (Rain)')
ylim([0,1])

figure(5)
boxplot(permute(errOn(4,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St5, 30 exp, 10 subexp, parameters box plot (Rain)')
ylim([0,1])

figure(6)
boxplot(permute(errOn(5,:,:),[3,2,1])), hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha$','v','h','$K_r$','$\gamma$','$K_p$','$\gamma_f$','$K_f$'})
xlabel('Parameter','Interpreter','Latex')
ylabel('${\epsilon}$','rotation',0,'Interpreter','Latex')
title('On-Line OED, St5, 30 exp, 20 subexp, parameters box plot (Rain)')
ylim([0,1])

% Test if the PE result converges in all the sub-exp
figure(7)
par_trajectory = [];
sub_exp_ind = 5;
for i = 1:numLoops(sub_exp_ind) 
    par = SubExp{sub_exp_ind,1}.pe_results(1,i);
    par = par{1,1}.fit.thetabest;
    par_trajectory(i,:) = par;
end
plot(par_trajectory)

%% Rain190123 Plot time vs parameter estimates
close all
cb = cbrewer('qual', 'Dark2', 30);
for i = 1:length(numLoops)
    for j = 1:numExperiments
        figure(i)
        h = plot(SubExp{i,j}.pe_results{1,1}.nlpsol.time, SubExp{i,j}...
            .pe_results{1,1}.nlpsol.bestit);
        set(h, {'color'}, {cb(j,:);  cb(j,:);  cb(j,:); cb(j,:);...
            cb(j,:); cb(j,:); cb(j,:); cb(j,:);})
        hold on
        plot([0, max(SubExp{i,j}.pe_results{1,1}.nlpsol.time)],[true_par',true_par'],...
            '--k','LineWidth',1.5);
        title(strcat(int2str(numLoops(i)), ' SubExp'))
        xlabel('Time')
    end
end
%%
% g = plot(pe_results{1,1}.nlpsol.time,pe_results{1,1}.nlpsol.bestit);
% hold on
% h = plot([0, max(pe_results{1,1}.nlpsol.time)],[true_par',true_par'],'--',...
%     'LineWidth',1.5);
% set(g, {'color'},{cb(1,:); cb(2,:); cb(3,:); cb(4,:); cb(5,:); cb(6,:);...
%     cb(7,:); cb(8,:)})
% set(h, {'color'},{cb(1,:); cb(2,:); cb(3,:); cb(4,:); cb(5,:); cb(6,:);...
%     cb(7,:); cb(8,:)})
% legend('alpha1','Vm1','h1','Km1','d1','alpha2','d2','Kf')
%% Rain190125 single parameter plot
cb = cbrewer('qual', 'Dark2', 30);
index = 1;
par_name = {'alpha1', 'Vm1', 'h1', 'Km1', 'd1', 'alpha2', 'd2', 'Kf'};
for i = 1:length(numLoops)
    for k = 1:length(true_par)
        figure
        for j = 1:numExperiments
            h = plot(SubExp{i,j}.pe_results{1,1}.nlpsol.time, SubExp{i,j}...
                .pe_results{1,1}.nlpsol.bestit(:,k));
            set(h, {'color'}, {cb(j,:)})
            hold on
            plot([0, max(SubExp{i,j}.pe_results{1,1}.nlpsol.time)],...
                [true_par(k)',true_par(k)'],...
                '--k','LineWidth',1.5);
            title(strcat(int2str(numLoops(i)), ' SubExp,', {' '}, par_name{k}))
            xlabel('Time')
        end
    end
end
%% Rain190218 Residual plot of the no noise experiment
% Rain_Step40_SubExp4_Exp50.mat is the experiment without noise
load('Rain_Step40_SubExp4_Exp50.mat')
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
[colormap] = cbrewer('qual','Set1',3);
subplot(2,1,1)
bar(best_global_theta'-true_par, 'FaceColor', colormap(2,:));
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha_1$','$Vm_1$','$h_1$','$Km_1$','$d_1$','$\alpha_2$',...
    '$d_2$','$Kf$'})
ylabel('Residual','Interpreter','Latex')
title('No noise','Interpreter','Latex')

% Rain_Step40_SubExp4_Exp30.mat is the experiment with noise of std = 0.05;
subplot(2,1,2)
load('Rain_Step40_SubExp4_Exp30.mat')
bar(best_global_theta'-true_par, 'FaceColor', colormap(2,:));
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(gca,'xticklabel',{'$\alpha_1$','$Vm_1$','$h_1$','$Km_1$','$d_1$','$\alpha_2$',...
    '$d_2$','$Kf$'})
ylabel('Residual','Interpreter','Latex')
title('Noise with std dev $0.05$','Interpreter','Latex')

% The parameter residuals of the no noise experiment are almost zero
% compare the noisy one.

