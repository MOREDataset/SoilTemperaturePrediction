%% MERGE SOIL DATA with MERRA-2 data
clc
clear
%original soil dataset
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik.csv");
opts = detectImportOptions(path1);
Air2UGalltimeseriesToolik = readtimetable(path1,opts);
Air2UGalltimeseriesToolik = sortrows(Air2UGalltimeseriesToolik);
dt = days(1);
Air2UGalltimeseriesToolik = retime(Air2UGalltimeseriesToolik,'regular','fillwithmissing','TimeStep',dt);

%MERRA-2 part1
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\ICOP2023',"MERRA2_dataprocessed.csv");
opts = detectImportOptions(path1);
MERRA2part1 = readtimetable(path1,opts);
MERRA2part1 = sortrows(MERRA2part1);
MERRA2part1 = retime(MERRA2part1,'regular','fillwithmissing','TimeStep',dt);

%MERRA-2 part2
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\ICOP2023',"MERRA2part2_dataprocessed.csv");
opts = detectImportOptions(path1);
MERRA2part2 = readtimetable(path1,opts);
MERRA2part2 = sortrows(MERRA2part2);
MERRA2part2 = retime(MERRA2part2,'regular','fillwithmissing','TimeStep',dt);
MERRA2part2.ALBEDO = [];

%MERRA-2 part3
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\ICOP2023',"MERRA2part3_dataprocessed.csv");
opts = detectImportOptions(path1);
MERRA2part3 = readtimetable(path1,opts);
MERRA2part3 = sortrows(MERRA2part3);
MERRA2part3 = retime(MERRA2part3,'regular','fillwithmissing','TimeStep',dt);

% Synchronize timetables
newTimetable = synchronize(Air2UGalltimeseriesToolik,MERRA2part1,MERRA2part2, MERRA2part3,...
"intersection","fillwithmissing");
newTimetable = rmmissing(timetable2table(newTimetable));
% newTimetable.GHTSKIN = [];
newTimetable.PRECTOT = [];
newTimetable.LAI = [];
newTimetable.GWETPROF = [];
newTimetable.GWETROOT = [];
newTimetable.GWETTOP = [];
% newTimetable.Tair = [];
% newTimetable.T0 = [];
newTimetable.T2M = newTimetable.T2M - 273.15;
newTimetable.TS = newTimetable.TS - 273.15;
newTimetable.TLML = newTimetable.TLML - 273.15;
newTimetable.TSH = newTimetable.TSH - 273.15;
% writetable(newTimetable,fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv"))
% % 
% % %original soil dataset with additional features from MERRA-2 dataset
% % path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv");
% % opts = detectImportOptions(path1);
% % tbl = readtimetable(path1,opts);
% % %original MERRA-2 part2 dataset
% % path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv");
% % opts = detectImportOptions(path1);
% % tbl = readtimetable(path1,opts);

dt = days(1);
tbl = retime(table2timetable(newTimetable),'regular','fillwithmissing','TimeStep',dt);
ftrs = {'TLML','SWLAND','SNODP','SWGDN','T2M'};
for i = 1:length(ftrs)
    fprintf('%s: %.2f & %.2f\n',ftrs{i},min(tbl.(ftrs{i})),max(tbl.(ftrs{i})))
end
% resols = [1,7,14,30];
% for J = 1:length(resols)
%     res = resols(J);
%     tbl = movingMeanStride(tbl, res);
%     RMSE = sqrt(mean((tbl.Tair-tbl.T2M).^2,"omitmissing"));
%     fprintf('res%d, AllAtOnce, RMSE=%.3f\n', res, RMSE)
%     %RESAMPLE    
%     seasons = {'Spring','Summer','Autumn','Winter'};
%     for S = 1:length(seasons)
%         season = seasons{S};
%         %generate the season specific table of data
%         TBL = season_splitting(tbl,season);
%         RMSE = sqrt(mean((TBL.Tair-TBL.T2M).^2,"omitmissing"));
%         fprintf('res%d, %s, RMSE=%.3f\n', res, season, RMSE)
%     end
% end
%% Season-based corr
clear
clc
res = 7;
num_top = 6;
trgts = {'T0','T16','T31','T46','T76','T97'};
seasons = {'Spring','Summer','Autumn','Winter'};
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv");
opts = detectImportOptions(path1);
tbl = readtimetable(path1,opts);
all_RemoteSensing_features = tbl.Properties.VariableNames(length(tbl.Properties.VariableNames)-15+1:end);
AVG_corr = table([strings(1,1),zeros(1,length(all_RemoteSensing_features))]);
AVG_corr = splitvars(AVG_corr, 'Var1');
AVG_corr.Properties.VariableNames = ['Season', all_RemoteSensing_features(:)'];
for S = 1:length(seasons)
    season = seasons{S};
    fprintf('%s:\n', season)
    for J = 1:length(trgts)
        depths = ["T0","T8","T16","T23","T31","T38","T46","T61","T76","T97","Tair"];
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv");
        opts = detectImportOptions(path1);
        tbl = readtimetable(path1,opts);
        dt = days(1);
        tbl = retime(tbl,'regular','fillwithmissing','TimeStep',dt);
        %RESAMPLE
        tbl = movingMeanStride(tbl, res);
        trgt = trgts{J};
        depths(depths == trgt) = [];
        tbl = removevars(tbl, depths);
        %generate the season specific table of data
        TBL = season_splitting(tbl,season);
        for K = 2:length(TBL.Properties.VariableNames)
            corrResults.(season).(trgt).(TBL.Properties.VariableNames{K}) = corr(TBL.(TBL.Properties.VariableNames{K}),TBL.(trgt));
        end
        [temp,temp_i] = maxk(abs(cellfun(@(x)(corrResults.(season).(trgt).(x)),fieldnames(corrResults.(season).(trgt)))),num_top);
        fprintf('- %s, top %d correlated features:', trgt, num_top)
        for ii = 1:num_top
            fprintf('%s,', string(TBL.Properties.VariableNames(temp_i(ii)+1)))
        end
        fprintf('\n')
    end
    tmp = struct2table([corrResults.(season).T0, corrResults.(season).T16, corrResults.(season).T31, corrResults.(season).T46, corrResults.(season).T76, corrResults.(season).T97]);
    AVG_corr.Season{S} = season;
    AVG_corr(S,2:end) = mean(tmp);
end
AVG_corr.Season = categorical(AVG_corr.Season);
tmp = 'AVG_corr.';
for I = 2:length(AVG_corr.Properties.VariableNames(:))
    if I == 2 
        tmpp = strcat(tmp,AVG_corr.Properties.VariableNames(I));
    else
        tmpp = strcat(tmpp,AVG_corr.Properties.VariableNames(I));
    end
    if I+1 <= length(AVG_corr.Properties.VariableNames(:))
        tmpp = strcat(tmpp, ',', tmp);
    end
    AVG_corr.(AVG_corr.Properties.VariableNames{I}) = double(AVG_corr.(AVG_corr.Properties.VariableNames{I}));
end

size1 = 30;
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout("flow",'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile, hold all
eval(sprintf('b = bar(AVG_corr.Season,[%s],''grouped'',''FaceColor'',''flat'');',string(tmpp{:})))
ylabel('Correlation coefficient','interpreter','latex','FontSize',size1);
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on')
yline(0.5,"LineStyle",'--',"LineWidth",2,'Label',"Threshold = 0.5",'Interpreter','latex',LabelHorizontalAlignment='right', LabelVerticalAlignment='middle',fontsize=size1-10)
yline(-0.5,"LineStyle",'--',"LineWidth",2,'Label',"Threshold = -0.5",'Interpreter','latex',LabelHorizontalAlignment='right', LabelVerticalAlignment='middle',fontsize=size1-10)
legend(string(AVG_corr.Properties.VariableNames(2:end)),'Location','northoutside','Interpreter','latex','Orientation','horizontal',FontSize=size1-5,NumColumns=5)
title(legend,'Remote sensing features','FontSize',size1,'FontWeight','bold','Interpreter','latex');
for k = 1:length(AVG_corr.Properties.VariableNames(2:end))
    b(k).CData = k;
end
cmap = colormap(gcf,sky);

% for k = 1:length(AVG_corr.Properties.VariableNames(2:end))
%     xtips1 = b(k).XEndPoints;
%     if b(k).YData>0
%         ytips1 = b(k).YEndPoints+.2;
%     else
%         ytips1 = b(k).YEndPoints-.2;
%     end
%     labels1 = string(round(b(k).YData,2));
%     text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
%         'VerticalAlignment','middle',FontSize=size1-10,Rotation=90)
% end

exportgraphics(gcf,strcat('CORRseasons_average6depths_ToolikLake.pdf'),"ContentType","image")


%% Season-based splitting for model learning

% Generate combined features and lags
clear
clc
lengths = [7*12];%[1,4,7,10,14,21,28,35,42,49,49,56,63,70]
resolution = 30;%[1,7,14,30];
res = strcat(num2str(resolution(1)),'D');
seasons = {'Spring','Summer','Autumn','Winter'};
trgts = {'T0','T16','T31','T46','T76','T97'};
for K = 1:length(lengths)
    LEN = lengths(K);
    fprintf('**************LEN=%d**********\n',LEN)
    for S = 1:length(seasons)
        season = seasons{S};
        for J = 1:length(trgts)
            trgt = trgts{J};
            path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv");
            opts = detectImportOptions(path1);
            tbl = readtimetable(path1,opts);
            dt = days(1);
            tbl = retime(tbl,'regular','fillwithmissing','TimeStep',dt);
            %RESAMPLE
            tbl = movingMeanStride(tbl, str2num(cell2mat(regexp(res,'\d*','match'))));
            feature_names = ["T0","T8","T16","T23","T31","T38","T46","T61","T76","T97"];
            feature_names(feature_names == trgt) = [];
            tbl = removevars(tbl, feature_names);
            names = tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,trgt));

            %generate the season specific table of data
            TBL = season_splitting(tbl,season);

            fprintf('%s, %s target, Length %d of features:', season, trgt, LEN)
            tbl2 = TBL(:,trgt);
            for f = 1:length(names)
                fprintf('%s, ', names{f})
                eval(sprintf('lags_%s = lagmatrix(TBL{:,names{%d}},0:LEN-1);',(names{f}),f))
                eval(sprintf('tbl2.%s = lags_%s;',(names{f}),(names{f})))
                tbl2 = splitvars(tbl2, (names{f}));
            end
            tbl2 = rmmissing(tbl2);
            if resolution==1
                RES= '';
            else
                RES = strcat('_res',res);
            end
            writetimetable(tbl2,fullfile(strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\','Air2UG_alltimeseries_Toolik_MERRA2data_',trgt,RES,'_SeqLen', num2str(LEN), '_',season,'.csv')))
            fprintf('\n')
            clear opts dt path1
        end
    end
end

%% * Compute Feature importance of LightGBDT
clc
size1=20;
path1 = "C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\ICOP2023\";
LEN = [14,21,28];
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,3,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
for i = 1:3
    nexttile,
    file_name = strcat('FeatureImportance_',num2str(LEN(i)),'_T76.csv');
    opts = detectImportOptions(fullfile(path1,file_name));% Import the data
    eval(sprintf('Importance_%s = readtable(fullfile(path1,file_name), opts);',num2str(LEN(i))))
    eval(sprintf('feature_names = categorical(Importance_%s.base_feature);',num2str(LEN(i))))
    eval(sprintf('feature_names = reordercats(feature_names, Importance_%s.base_feature);',num2str(LEN(i))))
    eval(sprintf('Importance_%s.base_feature = categorical(Importance_%s.base_feature);',num2str(LEN(i)),num2str(LEN(i))))
    eval(sprintf('p = bar(feature_names(4:end), Importance_%s.varimp(4:end));',num2str(LEN(i))))
    xtips1 = p(1).XEndPoints;
    ytips1 = p(1).YEndPoints+.75;
    labels1 = string(round(p(1).YData,2));
    title(sprintf('Using features of %s in length', num2str(LEN(i))))
    text(xtips1,ytips1,labels1,'HorizontalAlignment', 'Center', 'VerticalAlignment', 'middle',Interpreter='latex', FontSize=size1-7)
    set(gca,'FontSize',size1,'TickLabelInterpreter','latex','Ylim',[0 45])
    clear opts
end
fprintf('Most important features >5 in importance:\n')
Importance_21.base_feature(Importance_21.varimp>5)
%% *read data
clear
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv");
opts = detectImportOptions(path1);
tbl = readtimetable(path1,opts);
dt = days(1);
tbl = retime(tbl,'regular','fillwithmissing','TimeStep',dt);
tbl = removevars(tbl, ["T0", "T8","T16","T23","T31","T38","T46","T61"]);
tbl.T2M = tbl.T2M - 273.15;
tbl.TS = tbl.TS - 273.15;
tbl.TLML = tbl.TLML - 273.15;
tbl.TSH = tbl.TSH - 273.15;
names = tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T76'));
clear opts dt path1

%*correlation with moving window T76
clc 
size1 = 30;
cmap = get(0, 'defaultaxescolororder');

%remove T97 from table
tbl = tbl(:,tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T97')));
hrz = 1;
LEN = 21;%7:7:(7*8);
fprintf('Single-horizon correlation: \n')

names = tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T76'));
fprintf('Length %d of features: ', LEN(end))
for f = 1:length(names)
    fprintf('%s, ', names{f})
    eval(sprintf('lags_%s = lagmatrix(tbl{:,names{%d}},0:LEN(end)-1);',(names{f}),f))
    tbl2 = tbl(:,'T76');
    eval(sprintf('tbl2.%s = lags_%s;',(names{f}),(names{f})))
    tbl2 = splitvars(tbl2, (names{f}));
    tbl2 = rmmissing(tbl2);
    X = tbl2{:,(tbl2.Properties.VariableNames(2:end))};
    Y = tbl2{:,'T76'};
    % Calculate the Spearman correlation
    eval(sprintf('[CORR.%s,pvalue] = corr(X, Y, ''Type'', ''Spearman'');',names{f}))
end


%% *compute mean corr and plot T76
size1=30;
cmap = get(0, 'defaultaxescolororder');

fn = fieldnames(CORR);%get the name of features
corrTBL = table(categorical(strings(numel(fn),1)), zeros(numel(fn),1),'VariableNames',{'Feature','MCC'});%mean corr. coef (MCC)%initialize the table: feature, corr coef
for k=1:numel(fn)
    corrTBL(k,:) = table(categorical(fn(k)),mean(CORR.(fn{k})));
end
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout("flow",'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile, hold all
% %add legend workaround:
% barh(corrTBL.Feature(1),[nan nan]) %legend entries
% legend({'Remote sensing', 'in-situ measurements'},'Location','northoutside','Orientation','horizontal');

%then plot actual data
b=barh(corrTBL.Feature, corrTBL.("MCC"),HandleVisibility='off');
hold on,
xline(0.5,'--k','Threshold=0.5','Interpreter','latex',FontSize=size1-10,LabelHorizontalAlignment='center', LabelVerticalAlignment='bottom',HandleVisibility='off')
xline(-0.5,'--k','Threshold=-0.5','Interpreter','latex',FontSize=size1-10,LabelHorizontalAlignment='center', LabelVerticalAlignment='bottom',HandleVisibility='off')

xlabel('Spearman''s correlation coefficient','interpreter','latex','FontSize',size1);
% ylabel('Features','interpreter','latex','FontSize',size1);

set(gca,'FontSize',size1,'TickLabelInterpreter','latex','Xlim',[-1 1])
set(gca,'XGrid','on','XMinorGrid','on')
box("on")

xtips1 = b(1).YEndPoints;
for i=1:length(xtips1)
    if xtips1(i)>0
        xtips1(i)=xtips1(i)-0.06;
    else
        xtips1(i)=xtips1(i)+0.005;
    end
end
ytips1 = b(1).XEndPoints;
labels1 = string(round(b(1).YEndPoints,2));
text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
    'VerticalAlignment','middle', Interpreter='latex', FontSize=size1-11,Color=[1 1 1])
b.FaceColor = 'flat';
b.CData(1,:) = cmap(2,:);
for i=2:length(b.CData)
    b.CData(i,:) = cmap(1,:);
end
% legend1 =legend(gca,'show');
% set(legend1,'Orientation','horizontal', 'Interpreter','latex','FontSize',size1-7,'Location','northoutside')
% exportgraphics(t,strcat('CORR_T76vsAllfeatures.pdf'),'ContentType','vector')
% fprintf('Done!\n')
% exportgraphics(t,strcat('CORR_T76vsAllfeatures.pdf'),'ContentType','vector')
% close(gcf)

%% *read data, compute corr, mean corr and plot T97
clear
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv");
opts = detectImportOptions(path1);
tbl = readtimetable(path1,opts);
dt = days(1);
tbl = retime(tbl,'regular','fillwithmissing','TimeStep',dt);
tbl = removevars(tbl, ["T0", "T8","T16","T23","T31","T38","T46","T61"]);
tbl.T2M = tbl.T2M - 273.15;
tbl.TS = tbl.TS - 273.15;
tbl.TLML = tbl.TLML - 273.15;
tbl.TSH = tbl.TSH - 273.15;
clear opts dt path1

tbl = tbl(:,tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T76')));
hrz = 1;
LEN = 7:7:(90);
fprintf('Single-horizon correlation: \n')

names = tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T97'));
for f = 1:length(names)
    fprintf('Feature: %s with %d day lags', names{f}, LEN(end))
    eval(sprintf('lags_%s = lagmatrix(tbl{:,names{%d}},0:LEN(end)-1);',(names{f}),f))
    tbl2 = tbl(:,'T97');
    eval(sprintf('tbl2.%s = lags_%s;',(names{f}),(names{f})))
    tbl2 = splitvars(tbl2, (names{f}));
    tbl2 = rmmissing(tbl2);
    X = tbl2{:,(tbl2.Properties.VariableNames(2:end))};
    Y = tbl2{:,'T97'};
    % Calculate the Spearman correlation
    eval(sprintf('[CORR.%s,pvalue] = corr(X, Y, ''Type'', ''Spearman'');',names{f}))
end

size1=30;
cmap = get(0, 'defaultaxescolororder');

fn = fieldnames(CORR);%get the name of features
corrTBL = table(categorical(strings(numel(fn),1)), zeros(numel(fn),1),'VariableNames',{'Feature','MCC'});%mean corr. coef (MCC)%initialize the table: feature, corr coef
for k=1:numel(fn)
    corrTBL(k,:) = table(categorical(fn(k)),mean(CORR.(fn{k})));
end
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout("flow",'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile, hold all
% %add legend workaround:
% barh(corrTBL.Feature(1),[nan nan]) %legend entries
% legend({'Remote sensing', 'in-situ measurements'},'Location','northoutside','Orientation','horizontal');

%then plot actual data
b=barh(corrTBL.Feature, corrTBL.("MCC"),HandleVisibility='off');
hold on,
xline(0.5,'--k','Threshold=0.5','Interpreter','latex',FontSize=size1-10,LabelHorizontalAlignment='center', LabelVerticalAlignment='bottom',HandleVisibility='off')
xline(-0.5,'--k','Threshold=-0.5','Interpreter','latex',FontSize=size1-10,LabelHorizontalAlignment='center', LabelVerticalAlignment='bottom',HandleVisibility='off')

xlabel('Spearman''s correlation coefficient','interpreter','latex','FontSize',size1);
% ylabel('Features','interpreter','latex','FontSize',size1);

set(gca,'FontSize',size1,'TickLabelInterpreter','latex','Xlim',[-1 1])
set(gca,'XGrid','on','XMinorGrid','on')
box("on")

xtips1 = b(1).YEndPoints;
for i=1:length(xtips1)
    if xtips1(i)>0
        xtips1(i)=xtips1(i)-0.06;
    else
        xtips1(i)=xtips1(i)+0.005;
    end
end
ytips1 = b(1).XEndPoints;
labels1 = string(round(b(1).YEndPoints,2));
text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
    'VerticalAlignment','middle', Interpreter='latex', FontSize=size1-11,Color=[1 1 1])
b.FaceColor = 'flat';
b.CData(1,:) = cmap(2,:);
for i=2:length(b.CData)
    b.CData(i,:) = cmap(1,:);
end
% legend1 =legend(gca,'show');
% set(legend1,'Orientation','horizontal', 'Interpreter','latex','FontSize',size1-7,'Location','northoutside')
exportgraphics(t,strcat('CORR_T97vsAllfeatures.pdf'),'ContentType','vector')
fprintf('Done!\n')
exportgraphics(t,strcat('CORR_T97vsAllfeatures.pdf'),'ContentType','vector')
close(gcf)
%% correlation averaged over all values!
clc
%plot 
size1 = 30;
cmap = get(0, 'defaultaxescolororder');
cmap(8:13,:) = [0.721 0.216 0.671
0.933 0.498 0.164
0.051 0.549 0.459
0.886 0.290 0.482
0.157 0.706 0.514%
0.741 0.882 0.416];

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout("flow",'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile,
hold on,
numbest = 15;

%remove T97 from table
tbl = tbl(:,tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T97')));
tbl = rmmissing(tbl);
LEN = 1;
hrz = 1;
correlations = 0;
pvalue = 0;
fprintf('Single-horizon correlation: \n')

names = tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T76'));
for f = 1:length(names)
    fprintf('Feature: %s\n', names{f})
    X = tbl{:,(names{f})};
    Y = tbl{:,'T76'};
    tmp = [X,Y];
    tmp = rmmissing(tmp);
    % Calculate the Spearman correlation
    [correlations(f),pvalue(f)] = corr(tmp(:,1), tmp(:,2), 'Type', 'Spearman');
end
bar(categorical(names),correlations)
xlabel('Features','interpreter','latex','FontSize',size1);
ylabel('Spearman coefficient','interpreter','latex','FontSize',size1);
set(gca,'FontSize',size1,'TickLabelInterpreter','latex','Ylim',[-1 1])
set(gca,'YGrid','on','YMinorGrid','on')
box("on")
legend1 =legend(gca,'show');
set(legend1,'Orientation','horizontal','Interpreter','latex','FontSize',size1,'Location','northoutside')

clc
[~,I] = sort(abs(correlations),"descend");
fprintf('feature, correlation:\n')
for i=1:length(I)
    fprintf('%s, %.2f\n', names{I(i)}, correlations(I(i))*100)
end
% exportgraphics(t,strcat('CORR_T76plot.pdf'),'ContentType','vector')
% fprintf('Done!\n')
% exportgraphics(t,strcat('CORR_T76plot.pdf'),'ContentType','vector')
% fprintf('Done!\n')
% exportgraphics(t,strcat('CORR_T76plot.pdf'),'ContentType','vector')
% close(gcf)

%% using relieff instead
tbl = tbl(:,tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T97')));
names = tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T76'));
[idx,weights] = relieff(tbl{:,names}, tbl{:,'T76'},20);
figure,
hold on,
for i = 1:length(idx)
    barh(categorical(names(idx(i))), weights(idx(i)))
end
xlabel('Predictor rank')
ylabel('Predictor importance weight')

%% * heatmap correlation of inputs
clear
size1 = 30;
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Toolik_MERRA2data.csv");
opts = detectImportOptions(path1);
tbl = readtimetable(path1,opts);
dt = days(1);
tbl = retime(tbl,'regular','fillwithmissing','TimeStep',dt);
tbl = removevars(tbl, ["Tair", "T8","T16","T23","T31","T38","T46","T61"]);
tbl.T2M = tbl.T2M - 273.15;
tbl.TS = tbl.TS - 273.15;
tbl.TLML = tbl.TLML - 273.15;
tbl.TSH = tbl.TSH - 273.15;
tbl = tbl(:,tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T97')));
tbl = rmmissing(tbl);
names = tbl.Properties.VariableNames(~strcmp(tbl.Properties.VariableNames,'T76'));

% Define the extreme color and center color
extremeColor = [ 0    0.4470    0.7410];  % Color for -1 and 1./255
centerColor = [1 1 1];  % Color for 0./255

% Number of steps in each half of the colormap
nSteps = 10;

% Create the colormap
firstHalf = [linspace(extremeColor(1), centerColor(1), nSteps)', ...
             linspace(extremeColor(2), centerColor(2), nSteps)', ...
             linspace(extremeColor(3), centerColor(3), nSteps)'];

secondHalf = [linspace(centerColor(1), extremeColor(1), nSteps)', ...
              linspace(centerColor(2), extremeColor(2), nSteps)', ...
              linspace(centerColor(3), extremeColor(3), nSteps)'];

customColormap = [firstHalf; secondHalf];

% Compute the correlation matrix
correlationMatrix = corr(tbl{:,names});
% Create the heatmap
figure('units','normalized','outerposition',[0 0 1 1])
heatmap(names, names, round(correlationMatrix,2), 'Colormap',customColormap, ...
    'ColorLimits', [-1, 1], ...
    'MissingDataColor', [0.8 0.8 0.8]); 

set(gca,'FontSize',size1-12)%,'TickLabelInterpreter','latex'
set(gca,'XTickLabelRotation',90)
set(gcf, 'DefaultAxesFontName', 'Times New Roman');

exportgraphics(gcf,strcat('CORR_inputsplot.pdf'))
fprintf('Done!\n')
exportgraphics(gcf,strcat('CORR_inputsplot.pdf'))
fprintf('Done!\n')
exportgraphics(gcf,strcat('CORR_inputsplot.pdf'))
close(gcf)

% Add title and axis labels
% title('Heatmap of Feature Correlations');
% xlabel('Features');
% ylabel('Features');
% 
% % Define a correlation threshold; features with correlation above this will be considered "highly correlated"
% threshold = 0.95;  % Change this value based on your specific needs
% 
% % Find pairs of highly correlated features
% highlyCorrelatedPairs = {};
% for i = 1:size(correlationMatrix, 1)
%     for j = i+1:size(correlationMatrix, 2)
%         if abs(correlationMatrix(i, j)) > threshold
%             highlyCorrelatedPairs = [highlyCorrelatedPairs; {names{i}, names{j}}];
%         end
%     end
% end
% 
% % Display the pairs of highly correlated features
% disp('Highly correlated pairs of features:');
% disp(highlyCorrelatedPairs);




%% functions


function TBL = season_splitting(tbl,season)
mo = month(tbl.Date);
switch season
    case 'Spring'
        isSpring = ismember(mo, [3,4,5]);
        TBL = tbl(isSpring,:);
    case 'Summer'
        isSummer = ismember(mo, [6,7,8]);
        TBL = tbl(isSummer,:);
    case 'Autumn'
        isAutumn = ismember(mo, [9,10,11]);
        TBL = tbl(isAutumn,:);
    case 'Winter'
        isWinter = ismember(mo, [12,1,2]);
        TBL = tbl(isWinter,:);
end
end

function AveragedData = movingMeanStride(timetableData, resolution)
%resolution is the number of days over which the mean will be computed
    % Check if the input is a timetable
    if ~istimetable(timetableData)
        error('Input must be a timetable.');
    end
    % Check if the table is regular
    if ~isregular(timetableData)
        error('Timetable must be regular.');
    end
    % Calculate number of rows in a week
    % Assuming the timetable has daily data
    rowsInWeek = resolution;
    
    % Initialize the resulting timetable
    % The number of weeks will be the total days minus 6 (since we need 7 days for the first week)
    numWeeks = height(timetableData) - (rowsInWeek-1);
    AveragedData = timetable(timetableData.Date(1:numWeeks) + caldays(floor(rowsInWeek/2))); % Starting timepoints for each week to be middle of that week
    AveragedData.Properties.DimensionNames{1} = 'Date';

    % Iterate over each column of the timetable (excluding the Time column)
    for colIdx = 1:width(timetableData)
        colName = timetableData.Properties.VariableNames{colIdx};
        colData = timetableData{:, colIdx};

        % Apply moving mean with unit stride
        weeklyData = movmean(colData, [0, rowsInWeek-1], 'omitnan', 'Endpoints', 'discard');

        % Store the result in the weeklyMean timetable
        AveragedData.(colName) = weeklyData;
    end
end