clear all; clc; close all;

n_GaP = 3.22;
NumTap = [12;16;20];
NumMirr = [0;2;4;6;8;12;15;20];
rawdata = zeros(length(NumMirr),4,length(NumTap)); % num_taper, num_mirror, lambda, maxQ, V, E_overlap

for j=1:length(NumTap)
    for i=1:length(NumMirr)
        load("NumTap"+num2str(NumTap(j))+"/SiVminus_GaPonDia_wy360_wz222_NumTap"+num2str(NumTap(j))...
            +"_SweepNumMirr0To20_"+num2str(i)+".mat");
        rawdata(i,1,j) = taper_len;
        rawdata(i,2,j) = num_mir;
        rawdata(i,3,j) = lambda_maxQ;
        rawdata(i,4,j) = maxQ;
    end
end

figure;
set(0,'defaultFigureColor','w')
for j=1:length(NumTap) 
    semilogy(rawdata(:,2,j),rawdata(:,4,j),'-o','LineWidth',2);
    hold on;
end
hold off
xlim([0 20]);
ylim([0 1e6]);
set(gca,'XTick',[0 5 10 15 20])
set(gca,'YTick',[0 1e2 1e3 1e4 1e5 1e6])
set(gca, 'FontSize', 16)
legend({'12 Taper holes','16 Taper holes','20 Taper holes'},...
    'Location','northwest')
xlabel('Number of mirror holes')
ylabel('Q factor')
