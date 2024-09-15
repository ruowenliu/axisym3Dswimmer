% Copy Right @ Ruowen Liu, May 2023
% Save all results here and plot comparison figure
% Updated in Sep, 2024

%% [nu, JE, Jdrag]

nu_JE_Jdrag_prolate = [
    0.500000, 6.316561, 1.160235
    0.550000, 4.885154, 1.105492
    0.600000, 3.859771, 1.062557
    0.650000, 3.099078, 1.028850
    0.700000, 2.517108, 1.002602
    0.750000, 2.059092, 0.982620
    0.800000, 1.688483, 0.968175
    0.850000, 1.379398, 0.959001
    0.900000, 1.111299, 0.955513
    0.950000, 0.861774, 0.959899
    1.000000, 0.500000, 1.000000
    ];

nu_JE_Jdrag_maxeffi = [
    0.499877, 10.589063, 1.175896
    0.549316, 7.000447, 1.107867
    0.600517, 5.789659, 1.065386
    0.649267, 4.143658, 1.027017
    0.700391, 3.320296, 1.000908
    0.748526, 2.477294, 0.979603
    0.798655, 1.931338, 0.965315
    0.849080, 1.513392, 0.956819
    0.899626, 1.167998, 0.954731
    0.949408, 0.876645, 0.959545
    1.000000, 0.500000, 1.000000
];

nu_JE_Jdrag_dragmin = [
    0.499933, 7.410853, 1.159464
    0.550014, 5.785206, 1.104138
    0.600049, 4.830780, 1.060375
    0.649966, 3.868583, 1.026001
    0.700033, 3.083887, 0.999202
    0.749922, 2.442337, 0.979120
    0.799981, 1.915446, 0.965006
    0.850024, 1.497402, 0.956677
    0.899949, 1.157881, 0.954314
    0.949980, 0.871795, 0.959540
    1.000000, 0.500000, 1.000000
];


%%% close figures
close all

%%% define colors

colorgrey = '#707070'; % grey: prolate
colorblack = '#000000'; % black: prolate
colorblue = '#0000a7'; % blue: maxeffi
colorgreen = '#77AC30'; % green: dragmin
colorred = '#A2142F'; % red
colororange = '#F28522'; % orange

defaultfontsize = 20;

choosemarkersize = 10;

%%% plot nu vs JE
figure(1)

p1_dragmin = plot(nu_JE_Jdrag_dragmin(:,1), nu_JE_Jdrag_dragmin(:,2),"diamond",'MarkerFaceColor',colorgreen,'color',colorgreen,'MarkerSize',choosemarkersize);
hold on
p1_maxeffi = plot(nu_JE_Jdrag_maxeffi(:,1), nu_JE_Jdrag_maxeffi(:,2),"pentagram",'MarkerFaceColor',colororange,'color',colororange,'MarkerSize',choosemarkersize+3);
p1_prolate = plot(nu_JE_Jdrag_prolate(:,1), nu_JE_Jdrag_prolate(:,2),"o",'MarkerFaceColor',colorblue,'color',colorblue,'MarkerSize',choosemarkersize-2);

legend([p1_maxeffi,p1_dragmin,p1_prolate],'Max Efficiency','Min Drag Force','Prolate Spheroid','Location','northeast','NumColumns',1,'Interpreter','latex');

set(gca,'FontSize',defaultfontsize,'TickLabelInterpreter','latex'); 

xlabel('Reduced Volume','FontSize',defaultfontsize,'Interpreter','latex')
x1 = 0.50-0.02; x2 = 1.0+0.02; xlim([x1,x2]); set(gca,'XTick',0.5:0.05:1);
tixX=get(gca,'XTick')'; set(gca,'XTickLabel',num2str(tixX,'%.2f'));

ylabel('Swimming Efficiency','FontSize',defaultfontsize,'Interpreter','latex')
y1 = 0; y2 = 11; ylim([y1,y2]); set(gca,'YTick',y1:1.0:y2);
tixY=get(gca,'YTick')'; set(gca,'YTickLabel',num2str(tixY,'%.1f'));

text(0.75,11,"(a)","Interpreter","latex","FontSize",defaultfontsize);

%%% plot nu vs Jdrag
figure(2)

p2_dragmin = plot(nu_JE_Jdrag_dragmin(:,1), nu_JE_Jdrag_dragmin(:,3),"diamond",'MarkerFaceColor',colorgreen,'color',colorgreen,'MarkerSize',choosemarkersize);
hold on
p2_maxeffi = plot(nu_JE_Jdrag_maxeffi(:,1), nu_JE_Jdrag_maxeffi(:,3),"pentagram",'MarkerFaceColor',colororange,'color',colororange,'MarkerSize',choosemarkersize+3);
p2_prolate = plot(nu_JE_Jdrag_prolate(:,1), nu_JE_Jdrag_prolate(:,3),"o",'MarkerFaceColor',colorblue,'color',colorblue,'MarkerSize',choosemarkersize-2);

p2_lineONE = plot([0,1.5],[1,1],'--','LineWidth',0.8,'Color',[0.2, 0.4, 0.2]);

set(gca,'FontSize',defaultfontsize,'TickLabelInterpreter','latex'); 

xlabel('Reduced Volume','FontSize',defaultfontsize,'Interpreter','latex')
x1 = 0.50-0.02; x2 = 1.0+0.02; xlim([x1,x2]); set(gca,'XTick',0.5:0.05:1);
tixX=get(gca,'XTick')'; set(gca,'XTickLabel',num2str(tixX,'%.2f'));

ylabel('Normalized Drag Force','FontSize',defaultfontsize,'Interpreter','latex')
y1 = 0.94; y2 = 1.19; ylim([y1,y2]); set(gca,'YTick',y1:0.02:1.18);
tixY=get(gca,'YTick')'; set(gca,'YTickLabel',num2str(tixY,'%.2f'));

text(0.75,1.19,"(b)","Interpreter","latex","FontSize",defaultfontsize);
legend([p2_maxeffi,p2_dragmin,p2_prolate],'Max Efficiency','Min Drag Force','Prolate Spheroid','Location','northeast','NumColumns',1,'Interpreter','latex');


%%% save figures

figure(1)
saveas(gcf, 'nu_vs_efficiency_noline','pdf')
saveas(gcf, 'nu_vs_efficiency_noline','epsc')

figure(2)
saveas(gcf, 'nu_vs_Jdrag_noline','pdf')
saveas(gcf, 'nu_vs_Jdrag_noline','epsc')

