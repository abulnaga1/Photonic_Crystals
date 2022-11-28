clear all; clc; close all;

dev_num = [6];
hx_list = zeros(size(dev_num));
hy_list = zeros(size(dev_num));

hx_design = 70;
hy_design = 136;
wy_design = 350;

num_taper = 16;
acav = 162;
amir = 180;

%Wy_list = [0];

for j = 1:length(dev_num)
    fitResult = [];
    fileList = 1:13;

    for i=1:length(fileList)-1
        load(['dev',num2str(dev_num(j)), '/', strcat(num2str(fileList(i)), 'to', num2str(fileList(i+1))),'.mat'])
            
        tmp_fitResult = zeros(length(hole_xcor),7);
        tmp_fitResult(:,1) = hole_xcor';
        tmp_fitResult(:,2) = hole_ycor';
        tmp_fitResult(:,3) = hole_hx';
        tmp_fitResult(:,4) = hole_hy';
        tmp_fitResult(:,6) = sigma';
        tmp_fitResult(:,7) = angle';
        tmp_fitResult = sortrows(tmp_fitResult);
        tmp = circshift(tmp_fitResult(:,1),length(hole_xcor)-1) - tmp_fitResult(:,1);
        tmp_fitResult(1:end-1,5) = tmp(1:end-1); %inter-hole spacing
        fitResult = cat(1,fitResult,tmp_fitResult);
    end

    hole_xcor_list = nonzeros(fitResult(:,1))';
    hole_ycor_list = nonzeros(fitResult(:,2))';
    hole_hx_list = nonzeros(fitResult(:,3))';
    hole_hy_list = nonzeros(fitResult(:,4))';
    hole_spacing_list = nonzeros(fitResult(:,5))';
    sigma_list = nonzeros(fitResult(:,6))';
    angle_list = nonzeros(fitResult(:,7))';


    one_side_holes = (num_taper-1)/2.0;
    idx_list = linspace(-1*one_side_holes,one_side_holes,num_taper);
    aper_a = (amir - acav)/(one_side_holes*one_side_holes - 0.25);
    aper_c = (amir - 4*acav*one_side_holes*one_side_holes)/(1-4*one_side_holes*one_side_holes);
    acav_list = zeros(size(idx_list));
    for i=1:length(idx_list)
        acav_list(i) = aper_a*idx_list(i)*idx_list(i) + aper_c;
    end
    for i=1:length(acav_list)-1
        acav_spacing(i) = acav_list(i)/2 +  acav_list(i+1)/2;
    end

    num_mirr = (length(hole_spacing_list)-num_taper+1)/2;
    amirr_spacing = ones(1,num_mirr)*amir;
    aper_spacing = [amirr_spacing acav_spacing amirr_spacing];


    figure;
    ind = [1:length(hole_spacing_list)]+0.5;
    scatter(ind,hole_spacing_list,100,'k')
    hold on
    plot(ind,aper_spacing,'LineWidth',2)
    xlim([1 length(hole_spacing_list)+1])
    %ylim([160 200])
    set(gca,'FontSize',16)
    xlabel('Hole #','FontSize',18)
    ylabel('Spacing (nm)','FontSize',18)
    legend('Measured from SEM','CAD design value')
    saveas(gcf,['dev',num2str(dev_num(j)),'/HolePos.png'])

    figure;
    plot(ind,abs(aper_spacing-hole_spacing_list),'LineWidth',2)
    set(gca,'FontSize',16)
    xlim([1 length(hole_spacing_list)+1])
    xlabel('Hole #','FontSize',18)
    ylabel('Difference (nm)','FontSize',18)
    saveas(gcf,['dev',num2str(dev_num(j)),'/HolePosDiff.png'])

    hx = hole_hx_list;
    hy = hole_hy_list;
    diff = sigma_list;

    figure;
    ind = 1:length(hx);
    errorbar(ind,hx,diff,'LineWidth',2)
    hold on;
    yline(hx_design,'k--','linewidth',2);
    legend('measured hx','desired hx', 'location', 'southeast')
    %legend boxoff;
    set(gca,'FontSize',16)
    xlim([1 length(hx)])
    ylim([min(hx_design, min(hx)) - 5, max(hx_design, max(hx))+5])
    xlabel('Hole #','FontSize',18)
    ylabel('hx (nm)','FontSize',18)
    title(['Average hx = ',num2str(mean(hx)), ' nm'],'FontSize',18)
    saveas(gcf,['dev',num2str(dev_num(j)),'/hx.png'])

    figure;
    ind = 1:length(hy);
    errorbar(ind,hy,diff,'LineWidth',2)
    hold on;
    yline(hy_design,'k--','linewidth',2);
    legend('measured hy','desired hy', 'location', 'southeast')
    %legend boxoff;
    set(gca,'FontSize',16)
    xlim([1 length(hy)])
    ylim([min(hy_design, min(hy)) - 5, max(hy_design, max(hy))+5])
    xlabel('Hole #','FontSize',18)
    ylabel('hy (nm)','FontSize',18)
    title(['Average hy = ',num2str(mean(hy)), ' nm'],'FontSize',18)
    saveas(gcf,['dev',num2str(dev_num(j)),'/hy.png'])
    
    close all;
    hx_list(j) = mean(hx);
    hy_list(j) = mean(hy);
end

%wy_list = [747;749;742].*200./480;

% figure;
% plot(dev_num,wy_list,'LineWidth',2,'marker','o');
% hold on
% yline(wy_design,'k-.','LineWidth',2);
% hold off
% set(gca,'FontSize',16)
% set(gca,'XTick',dev_num);
% ylim([min(wy_list)-10; wy_design+10])
% xlabel('Dev #','FontSize',18)
% ylabel('measured value of wy (nm)','FontSize',18)
% legend('measured wy','desired wy')
% saveas(gcf,'wy_vs_design.png')

figure;
plot(dev_num,hx_list,'LineWidth',2,'marker','o');
hold on
plot(dev_num,hx_design*ones(size(dev_num)),'k-.','LineWidth',2);
%yline(hx_design,'k--','linewidth',2);
hold off
set(gca,'FontSize',16)
%ylim([hx_design-10, max(hx_list)+5]);
set(gca,'XTick',dev_num);
xlabel('Dev #','FontSize',18)
ylabel('measured value of hx (nm)','FontSize',18)
legend('measured hx','desired hx','location','southeast')
saveas(gcf,'hx_vs_design.png')

figure;
plot(dev_num,hy_list,'LineWidth',2,'marker','o');
hold on
plot(dev_num,hy_design*ones(size(dev_num)),'k-.','LineWidth',2);
%yline(hy_design,'k--','linewidth',2);
hold off
set(gca,'FontSize',16)
%ylim([hx_design-10, max(hy_list)+5]);
set(gca,'XTick',dev_num);
xlabel('Dev #','FontSize',18)
ylabel('measured value of hy (nm)','FontSize',18)
legend('measured hy','desired hy')
saveas(gcf,'hy_vs_design.png')

