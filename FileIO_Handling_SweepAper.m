% Reading the results of full band structure calculations from MPB

clc; close all; clear all;

% Structure of output.dat: 
% symmetry, folderID, wy, wz, aspect, aper, ff, k_index, kx, ky, kz, kmag, 
% band 1, band 2

fid = fopen('output.dat','r');
output_result = textscan(fid,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f',...
    'Delimiter',',');
fclose(fid);

set(0,'defaultFigureColor','w')

neff = 2.6;
lambda_adjusted = 1221;
aper= output_result{6};
band1= output_result{13};
band2= output_result{14};

band1 = band1(1:2:end); 
band2 = band2(1:2:end);
aper = aper(1:2:end);

scatter(aper,band1*(3e5)./aper)
hold on
scatter(aper,band2*(3e5)./aper)
scatter(aper,(band1+band2)*(3e5)./(2*aper))
%scatter(aper,band3*2*(3e17)*pi./aper)
%scatter(aper,band4*2*(3e17)*pi./aper)
%plot(sort(aper),((3e5)/lambda_adjusted)*ones(size(aper)),'k--','LineWidth',2)
%plot([150 150],[2.5e2 6.5e2],'k-','LineWidth',2)
%plot([177 177],[2.5e2 6.5e2],'k-.','LineWidth',2)
%plot(aper,(324.45071e12)*ones(size(aper)),'k-.','LineWidth',2)
hold off
set(gca, 'FontSize', 16)
xlabel('aper (nm)','FontSize',22)
ylabel('Frequency (THz)','FontSize',22)
title('wy=350nm, wz=220nm','FontSize',22)
legend('1st band edge', '2nd band edge', 'mid-gap','designed cavity resonance','a_{cav} = 206nm', 'a_{mirr} = 253nm')
%axis([130 190 250 450])