clc,close all;
clear all;
warning off all;
noOfNodes=100;
rounds=2500;
XYarea=100;
%% method
% 0:    Hardware with compression
% 1:    Software no compression (Traditional)
% 2:    Software with compression
%% Simulate when hardware is used to implement functions
disp('Simulation when key functions are implemented in hardware');
method=0;
Compression = 1;
[ResultsPerRoundH,sensorsPerCHperRoundH,CHperRoundH, maxBlkSizePerRoundH,...
    minBlkSizePerRoundH,TotalExtraFileSizeinBytesH,LifeExpectancyH,...
    SensorsHardwareH]=simulateWSN(method,Compression,noOfNodes,XYarea,rounds);
disp('End simulation when key functions are implemented in hardware');

%% Simulate when Software is used to implement functions
disp('Simulation when everything is by software without compression');
method=1;
Compression=0;
[ResultsPerRoundS,sensorsPerCHperRoundS,CHperRoundS, maxBlkSizePerRoundS,...
    minBlkSizePerRoundS,TotalExtraFileSizeinBytesS,LifeExpectancyS,...
    SensorsSoftware,CompRatio]=simulateWSN(method,Compression,noOfNodes,XYarea,rounds);
disp('End simulation when everything is by software without compression');
%% Simulate when software with compression
disp('Simulation when software with compression');
method=2;
Compression = 1;
[ResultsPerRoundSWC,sensorsPerCHperRoundSWC,CHperRoundSWC,...
    maxBlkSizePerRoundSWC, minBlkSizePerRoundSWC,...
    TotalExtraFileSizeinBytesSWC,LifeExpectancySWC,...
    SensorsHardwareSWC]=simulateWSN(method,Compression,noOfNodes,XYarea,rounds);
disp('End simulation when software with compression');
%% Unify sizes
[r1,~]=size(ResultsPerRoundS);
[r2,~]=size(ResultsPerRoundH);
[r3,~]=size(ResultsPerRoundSWC);
r=min([r1,r2,r3]);
%% Check Energy Differences HW with traditional
change=(ResultsPerRoundH.TotalConsumedEnergy(1:r)-ResultsPerRoundS.TotalConsumedEnergy(1:r));
percentGain = 100*mean(abs(change))/mean(ResultsPerRoundS.TotalConsumedEnergy(1:r));
%% Average Complexity HW with traditional
change = mean(ResultsPerRoundS.TotalComplexity)-mean(ResultsPerRoundH.TotalComplexity);
avgComplexityGain = 100*mean(change./mean(ResultsPerRoundS.TotalComplexity));

%% Check Energy Differences Software with compression with traditional
change=(ResultsPerRoundSWC.TotalConsumedEnergy(1:r)-ResultsPerRoundS.TotalConsumedEnergy(1:r));
percentGainSWC = 100*mean(abs(change))/mean(ResultsPerRoundS.TotalConsumedEnergy(1:r));
%% Average Complexity HW with traditional
change = mean(ResultsPerRoundS.TotalComplexity)-mean(ResultsPerRoundSWC.TotalComplexity);
avgComplexityGainSWC = 100*mean(change./mean(ResultsPerRoundS.TotalComplexity));

% %% Memory Calculations
% memorySizeHaedware = zeros(noOfNodes,r);
% memorySizeSoftware = zeros(noOfNodes,r);
% for sen=1:noOfNodes
%     memorySizeHaedware(sen,:) = (SensorsHardwareH(sen).memoryPerRound(1:r));
%     memorySizeSoftware(sen,:) = (SensorsSoftware(sen).memoryPerRound(1:r));
% end
% 
% 
% averageMemoryPerRoundHardware = mean(memorySizeHaedware);
% averageMemoryPerRoundSoftware = mean(memorySizeSoftware);
% averageMemoryPerRoundSoftwarePlus = mean(memorySizeSoftware)+TotalExtraFileSizeinBytesH ;
% 
% averageMemoryPerRound = [averageMemoryPerRoundSoftware;averageMemoryPerRoundHardware];
% averageMemoryPerRoundPlus = [averageMemoryPerRoundSoftwarePlus;averageMemoryPerRoundHardware];
% 
% %% Average Memory Gain
% change = averageMemoryPerRoundSoftwarePlus - averageMemoryPerRoundHardware;
% avgMemoryGain = 100*mean(change./averageMemoryPerRoundSoftwarePlus);
%% Energy Plots
figure('Name','Total Energy Per Round');
plot(ResultsPerRoundH.TotalEnergy);
hold on;
plot(ResultsPerRoundS.TotalEnergy);
hold off;
grid on;
xlabel('Rounds');
ylabel('Energy in Joules');
legend('Hardware compression','Traditional Software');
title(['Total Energy Gain %= ', num2str(percentGain)]);

figure('Name','Total Complexity Per Round');
plot(ResultsPerRoundH.TotalComplexity);
hold on;
plot(ResultsPerRoundS.TotalComplexity);
hold off;
grid on;
xlabel('Rounds');
ylabel('Time Complexity in ms');
legend('Hardware','Software');
title(['Average Complexity Gain %= ', num2str(avgComplexityGain), ...
        '  Extra Size=', ...
        num2str(TotalExtraFileSizeinBytesH),'  Bytes']);

figure('Name','Alive Nodes Per Round');
plot(ResultsPerRoundH.AliveSensors);
hold on;
plot(ResultsPerRoundS.AliveSensors);
hold off;
grid on;
xlabel('Rounds');
ylabel('No of alive nodes');
legend(['Hardware:Life Expect= ',num2str(LifeExpectancyH)],['Software, Life Expec=', num2str(LifeExpectancyS)]);
title(['Alive Nodes in the network ']);

figure('Name','Energy Variance');
plot(ResultsPerRoundH.EnergyVariance);
hold on;
plot(ResultsPerRoundS.EnergyVariance);
hold off;
grid on;
xlabel('Rounds');
ylabel('Energy Variance \sigma^2');
legend('Hardware','Software');
title(['Energy Variance per round']);

SensorToPlot=7;
figure('Name','Energy for a single sensor');
plot(SensorsHardwareH(SensorToPlot).EpR);
hold on;
plot(SensorsSoftware(SensorToPlot).EpR);
hold off;
grid on;
xlabel('Rounds');
ylabel('Energy in  Joules');
ylim([0 0.5]);
legend('Hardware','Software');
title(['Energy per round for sensor ', num2str(SensorToPlot)]);


%% Energy Plots between software with compression and traditional
figure('Name','Software with Compression Total Energy Per Round');
plot(ResultsPerRoundSWC.TotalEnergy);
hold on;
plot(ResultsPerRoundS.TotalEnergy);
hold off;
grid on;
xlabel('Rounds');
ylabel('Energy in Joules');
legend('SW w/C','Software');
title(['SW w/C Total Energy Gain %= ', num2str(percentGainSWC)]);

figure('Name','Software with Compression Total Complexity Per Round');
plot(ResultsPerRoundSWC.TotalComplexity);
hold on;
plot(ResultsPerRoundS.TotalComplexity);
hold off;
grid on;
xlabel('Rounds');
ylabel('Time Complexity in ms');
legend('SW w/C','Software');
title(['SW w/C Average Complexity Gain %= ', num2str(avgComplexityGainSWC), ...
        '  Extra Size=', ...
        num2str(TotalExtraFileSizeinBytesSWC),'  Bytes']);

figure('Name','Software with Compression Alive Nodes Per Round');
plot(ResultsPerRoundSWC.AliveSensors);
hold on;
plot(ResultsPerRoundS.AliveSensors);
hold off;
grid on;
xlabel('Rounds');
ylabel('No of alive nodes');
legend(['SW w/C:Life Expect= ',num2str(LifeExpectancySWC)],['Software, Life Expec=', num2str(LifeExpectancyS)]);
title(['Software with Compression Alive Nodes in the network ']);

figure('Name','Software with Compression Energy Variance');
plot(ResultsPerRoundSWC.EnergyVariance);
hold on;
plot(ResultsPerRoundS.EnergyVariance);
hold off;
grid on;
xlabel('Rounds');
ylabel('Energy Variance \sigma^2');
legend('SW w/Compress','Software');
title(['  Software with Compression Energy Variance per round']);

SensorToPlot=7;
figure('Name','Software with Compression Energy for a single sensor');
plot(SensorsHardwareSWC(SensorToPlot).EpR);
hold on;
plot(SensorsSoftware(SensorToPlot).EpR);
hold off;
grid on;
xlabel('Rounds');
ylabel('Energy in  Joules');
ylim([0 0.5]);
legend('Software w/compression','Software');
title(['Software with Compression Energy per round for sensor ', num2str(SensorToPlot)]);

%% Compare three methods
figure('Name','Total Complexity Per Round for three methods');
plot(ResultsPerRoundH.TotalComplexity);
hold on;
plot(ResultsPerRoundSWC.TotalComplexity);
hold on;
plot(ResultsPerRoundS.TotalComplexity);
hold off;
grid on;
xlabel('Rounds');
ylabel('Time Complexity in ms');
legend('HW', 'SW w/C','Traditional');
title('Total Complexity Per Round for three methods');



%% Memory Plot
% SensorToPlot=7;
% figure('Name','Memory for a single sensor');
% plot(SensorsHardwareH(SensorToPlot).memoryPerRound);
% hold on;
% plot(SensorsSoftware(SensorToPlot).memoryPerRound);
% hold off;
% grid on;
% xlabel('Rounds');
% ylabel('Memory in  Bytes');
% legend('Hardware','Software');
% title(['Memory per round for sensor ', num2str(SensorToPlot)]);
% 
% figure('Name','Average Memory for the WSN per round');
% bar(averageMemoryPerRound');
% grid on;
% le = legend('Software','Hardware');
% set(le,'location','northwest');
% grid on;
% title('Average Memory Per Round Comparison');
% xlabel('Rounds');
% ylabel('Size in Bytes');
% 
% figure('Name','Average Memory for the WSN per round');
% bar((averageMemoryPerRoundPlus'));
% legend('Software','Hardware');
% grid on;
% title({
%     ['Average Memory Per Round Comparison Plus Program memory.']
%     ['Average Memory Gain = ', num2str(avgMemoryGain),'%']
%     });
% xlabel('Rounds');
% ylabel('Size in Bytes');
% 
% 
% %% Memory Plots
% maxMemroyPerNodeHardware = max(memorySizeHaedware');
% figure,bar(maxMemroyPerNodeHardware');
% grid on;
% title({
%     ['Hardware Maximum Memory Size Per Node,']
%     [' Max Memory for WSN = ', num2str(max(maxMemroyPerNodeHardware))]
%     });
% xlabel('Nodes');
% ylabel('Size in Bits');
% 
% maxMemroyPerNodeSoftware = max(memorySizeSoftware')+TotalExtraFileSizeinBytesS;
% figure,bar((maxMemroyPerNodeSoftware')+TotalExtraFileSizeinBytesS);
% grid on;
% title({
%     ['Software Maximum Memory Size Per Node,']
%     [' Max Memory for WSN = ', num2str(max(maxMemroyPerNodeSoftware))]
%     });
% xlabel('Nodes');
% ylabel('Size in Bits');
% 
% %% CH Plot
% 
% figure('Name','CH Per Round fo Sotware');
% bar(CHperRoundS);
% grid on;
% xlabel('Rounds');
% ylabel('No of Cluster Heads');
% title({['Min Block Size: ', num2str(minBlkSizePerRoundS(2))]
%     ['Max Block Size: ', num2str(maxBlkSizePerRoundS(2))]});
% % hold on;
% % bar(CHperRoundS);
% % hold off;
% % legend('Hardware','Software');
% % title({['Average CH per round:']
% %      ['for Hardware: ',num2str(mean(CHperRound))]
% %      ['for Software: ',num2str(mean(CHperRoundS))]});
% 
% figure('Name','CH Per Round Hardware');
% bar(CHperRoundH);
% grid on;
% xlabel('Rounds');
% ylabel('No of Cluster Heads');
% title({['Min Block Size: ', num2str(minBlkSizePerRoundH(2))]
%     ['Max Block Size: ', num2str(maxBlkSizePerRoundH(2))]});
% 
