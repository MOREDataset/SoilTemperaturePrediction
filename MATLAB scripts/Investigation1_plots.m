%% Plot X: investigate best resolution using best ML per season - Deadhorse
Winter_RMSE = [5.06	3.94	3.43	3.05	4.38
4.17	3.47	3.19	2.91	4.22
2.90	3.07	2.60	2.45	3.59
2.09	1.83	1.81	1.63	2.09
1.61	1.75	1.70	1.33	1.82
1.57	1.56	1.57	1.60	1.62
];

Spring_RMSE = [2.77	1.75	1.40	1.34	2.47
2.80	2.18	1.89	1.85	3.19
3.07	2.56	2.28	1.71	3.01
1.68	1.48	1.32	1.32	2.18
1.63	1.50	1.19	1.31	1.97
1.48	1.43	1.31	1.26	1.83
];

Summer_RMSE = [2.13	1.01	1.10	0.63	0.99
2.06	1.49	1.29	0.92	1.09
2.52	1.32	0.94	1.08	1.45
1.49	1.18	1.23	0.99	1.57
0.33	0.34	0.35	0.32	0.32
0.34	0.39	0.43	0.41	0.31
];

Fall_RMSE = [3.79	3.00	2.10	2.54	2.78
3.60	3.82	3.20	2.86	3.65
2.54	2.22	1.91	1.72	2.45
0.59	0.56	0.51	0.47	0.61
0.15	0.14	0.07	0.07	0.19
0.03	0.03	0.05	0.06	0.02
];

seasons = {'Daily','Weekly','Biweekly','Monthly'};
size1 = 35;

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile,
bar(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), Winter_RMSE(:,1:4),'EdgeColor','none',FaceAlpha= 0.6);
hold on 
p = plot(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), Winter_RMSE(:,5),"LineStyle",":","Marker","o", 'LineWidth',2);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}));
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on', 'XTickLabelRotation',90)%,'Ylim',[-1 1]
title('Winter','FontSize',size1,'Interpreter','latex')

leg2=legend(gca,p,'Baseline','Location','northeast');set(leg2,'Interpreter','latex','Orientation','horizontal',FontSize=size1-5);

nexttile,
bar(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), Spring_RMSE(:,1:4),'EdgeColor','none',FaceAlpha= 0.6);
hold on 
p = plot(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), Spring_RMSE(:,5),"LineStyle",":","Marker","o", 'LineWidth',2);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}));
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on', 'XTickLabelRotation',90)%,'Ylim',[-1 1]
title('Spring','FontSize',size1,'Interpreter','latex')

leg2=legend(gca,p,'Baseline','Location','northeast');set(leg2,'Interpreter','latex','Orientation','horizontal',FontSize=size1-5);


nexttile,
bar(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), Summer_RMSE(:,1:4),'EdgeColor','none',FaceAlpha= 0.6);
hold on 
p = plot(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), Summer_RMSE(:,5),"LineStyle",":","Marker","o", 'LineWidth',2);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}));
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on', 'XTickLabelRotation',90)%,'Ylim',[-1 1]
title('Summer','FontSize',size1,'Interpreter','latex')

leg2=legend(gca,p,'Baseline','Location','northeast');set(leg2,'Interpreter','latex','Orientation','horizontal',FontSize=size1-5);

nexttile,
bar(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), Fall_RMSE(:,1:4),'EdgeColor','none',FaceAlpha= 0.6);
hold on 
p = plot(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), Fall_RMSE(:,5),"LineStyle",":","Marker","o", 'LineWidth',2);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}));
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on', 'XTickLabelRotation',90)%,'Ylim',[-1 1]
title('Fall','FontSize',size1,'Interpreter','latex')

colororder("gem")

ylabel(t,'RMSE ($^{\circ}$C)','interpreter','latex','FontSize',size1);
xlabel(t,'Soil Temperatures','interpreter','latex','FontSize',size1);

lgd = legend({'$Daily$','$Weekly$','$Biweekly$','$Monthly$'}, 'Location','layout','Interpreter','latex','Orientation','horizontal',FontSize=size1-5,Parent=t);
title(lgd,'Resolution')
lgd.Layout.Tile = 'North'; % <----- relative to tiledlayout

ah1=axes('position',get(gca,'position'),'visible','off');
leg2=legend(gca,p,'Baseline','Location','northeast');set(leg2,'Interpreter','latex','Orientation','horizontal',FontSize=size1-5);

exportgraphics(t,'Investigation1_BestResolution_Deadhorse.pdf',"ContentType","image")
close("all")
%% Plot X: investigate best resolution using best ML per season - Toolik lake

Winter_RMSE = [5.06	3.94	3.43	3.05	4.38
4.17	3.47	3.19	2.91	4.22
2.90	3.07	2.60	2.45	3.59
2.09	1.83	1.81	1.63	2.09
1.61	1.75	1.70	1.33	1.82
1.57	1.56	1.57	1.60	1.62
];

Spring_RMSE = [2.77	1.75	1.40	1.34	2.47
2.80	2.18	1.89	1.85	3.19
3.07	2.56	2.28	1.71	3.01
1.68	1.48	1.32	1.32	2.18
1.63	1.50	1.19	1.31	1.97
1.48	1.43	1.31	1.26	1.83
];

Summer_RMSE = [2.13	1.01	1.10	0.63	0.99
2.06	1.49	1.29	0.92	1.09
2.52	1.32	0.94	1.08	1.45
1.49	1.18	1.23	0.99	1.57
0.33	0.34	0.35	0.32	0.32
0.34	0.39	0.43	0.41	0.31
];

Fall_RMSE = [3.79	3.00	2.10	2.54	2.78
3.60	3.82	3.20	2.86	3.65
2.54	2.22	1.91	1.72	2.45
0.59	0.56	0.51	0.47	0.61
0.15	0.14	0.07	0.07	0.19
0.03	0.03	0.05	0.06	0.02
];

seasons = {'Daily','Weekly','Biweekly','Monthly'};
size1 = 35;

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile,
bar(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), Winter_RMSE(:,1:4),'EdgeColor','none',FaceAlpha= 0.6);
hold on 
p = plot(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), Winter_RMSE(:,5),"LineStyle",":","Marker","o", 'LineWidth',2);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}));
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on', 'XTickLabelRotation',90)%,'Ylim',[-1 1]
title('Winter','FontSize',size1,'Interpreter','latex')

leg2=legend(gca,p,'Baseline','Location','northeast');set(leg2,'Interpreter','latex','Orientation','horizontal',FontSize=size1-5);

nexttile,
bar(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), Spring_RMSE(:,1:4),'EdgeColor','none',FaceAlpha= 0.6);
hold on 
p = plot(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), Spring_RMSE(:,5),"LineStyle",":","Marker","o", 'LineWidth',2);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}));
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on', 'XTickLabelRotation',90)%,'Ylim',[-1 1]
title('Spring','FontSize',size1,'Interpreter','latex')

leg2=legend(gca,p,'Baseline','Location','northeast');set(leg2,'Interpreter','latex','Orientation','horizontal',FontSize=size1-5);


nexttile,
bar(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), Summer_RMSE(:,1:4),'EdgeColor','none',FaceAlpha= 0.6);
hold on 
p = plot(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), Summer_RMSE(:,5),"LineStyle",":","Marker","o", 'LineWidth',2);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}));
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on', 'XTickLabelRotation',90)%,'Ylim',[-1 1]
title('Summer','FontSize',size1,'Interpreter','latex')

leg2=legend(gca,p,'Baseline','Location','northeast');set(leg2,'Interpreter','latex','Orientation','horizontal',FontSize=size1-5);

nexttile,
bar(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), Fall_RMSE(:,1:4),'EdgeColor','none',FaceAlpha= 0.6);
hold on 
p = plot(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), Fall_RMSE(:,5),"LineStyle",":","Marker","o", 'LineWidth',2);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}));
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on', 'XTickLabelRotation',90)%,'Ylim',[-1 1]
title('Fall','FontSize',size1,'Interpreter','latex')

colororder("gem")

ylabel(t,'RMSE ($^{\circ}$C)','interpreter','latex','FontSize',size1);
xlabel(t,'Soil Temperatures','interpreter','latex','FontSize',size1);

lgd = legend({'$Daily$','$Weekly$','$Biweekly$','$Monthly$'}, 'Location','layout','Interpreter','latex','Orientation','horizontal',FontSize=size1-5,Parent=t);
title(lgd,'Resolution')
lgd.Layout.Tile = 'North'; % <----- relative to tiledlayout

ah1=axes('position',get(gca,'position'),'visible','off');
leg2=legend(gca,p,'Baseline','Location','northeast');set(leg2,'Interpreter','latex','Orientation','horizontal',FontSize=size1-5);

exportgraphics(t,'Investigation1_BestResolution_Toolik.pdf',"ContentType","image")
close("all")
