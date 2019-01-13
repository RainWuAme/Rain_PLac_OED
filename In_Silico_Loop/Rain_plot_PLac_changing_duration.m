clear all
true_par = [0.0164186333380725 0.291556643109224 1.71763487775568 ...
    5.14394334860864 0.229999999999978 6.63776658557266...
    0.00575139649497780 0.0216999999961899];
cd('D:\rain\PLac_result\Data_change_duration')
Sampling_time_on_line = 5;
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
