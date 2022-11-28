clear all; clc; close all;

n_GaP = 3.22;

figure;
load('SiVminus_SuspendedGaP_wy360_wz222_NumTap12_NumMirr6_transmission_1.mat');
lambda = MEM.lambda;
plot(lambda*(1e9),Tnet,'LineWidth',2)
set(gcf,'color','w');
set(gca, 'FontSize', 16)
xlabel('Wavelength (nm)')
ylabel('Normalized transmission')
xlim([500 900])

figure;
load('SiVminus_GaPonDia_wy360_wz222_NumTap12_NumMirr6_transmission_1.mat');
lambda = MEM.lambda;
plot(lambda*(1e9),Tnet,'LineWidth',2)
set(gcf,'color','w');
set(gca, 'FontSize', 16)
xlabel('Wavelength (nm)')
ylabel('Normalized transmission')
xlim([500 900])
