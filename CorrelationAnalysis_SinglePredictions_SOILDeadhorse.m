%% MERGE Deadhorse SOIL DATA with MERRA-2 data
clc
clear
%original soil dataset
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Deadhorse.csv");
opts = detectImportOptions(path1);
Air2UGalltimeseriesDeadhorse = readtimetable(path1,opts);
Air2UGalltimeseriesDeadhorse = sortrows(Air2UGalltimeseriesDeadhorse);
dt = days(1);
Air2UGalltimeseriesDeadhorse = retime(Air2UGalltimeseriesDeadhorse,'regular','fillwithmissing','TimeStep',dt);

%MERRA-2
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\ICOP2023',"MERRA2_Deadhorse_dataprocessed.csv");
opts = detectImportOptions(path1);
MERRA2 = readtimetable(path1,opts);
MERRA2 = sortrows(MERRA2);
MERRA2 = retime(MERRA2,'regular','fillwithmissing','TimeStep',dt);

% Synchronize timetables
newTimetable = synchronize(Air2UGalltimeseriesDeadhorse,MERRA2,...
"intersection","fillwithmissing");
newTimetable = rmmissing(timetable2table(newTimetable));
newTimetable.GWETPROF = [];
newTimetable.GWETROOT = [];
newTimetable.GWETTOP = [];
% newTimetable.GHTSKIN = [];
% newTimetable.Tair = [];
% newTimetable.T0 = [];
newTimetable.T2M = newTimetable.T2M - 273.15;
newTimetable.TS = newTimetable.TS - 273.15;
newTimetable.TLML = newTimetable.TLML - 273.15;
newTimetable.TSH = newTimetable.TSH - 273.15;
newTimetable = movevars(newTimetable, "SNODP", "Before", "GHTSKIN");
newTimetable = movevars(newTimetable, "SWGDN", "Before", "GHTSKIN");
newTimetable = movevars(newTimetable, "LWGAB", "Before", "GHTSKIN");
newTimetable = movevars(newTimetable, "T2M", "Before", "GHTSKIN");
newTimetable = movevars(newTimetable, "SWLAND", "Before", "GHTSKIN");

start_date = datetime(1900,1,1);
end_date = datetime(2000,05,30);
S = timerange(start_date,end_date,'closed');
newTimetable = table2timetable(newTimetable);
newTimetable = newTimetable(S,:);
% writetimetable(newTimetable,fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Deadhorse_MERRA2data.csv"))
newTimetable = timetable2table(newTimetable);
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
trgts = {'T0','T_12','T_22','T_32','T_62','T_72'};
seasons = {'Spring','Summer','Autumn','Winter'};
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Deadhorse_MERRA2data.csv");
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
        depths = ["T0","T_02","T_07","T_12","T_22","T_32","T_42","T_62","T_67","T_72","Tair"];
        path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Deadhorse_MERRA2data.csv");
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
            corrResults.(season).(trgt).(TBL.Properties.VariableNames{K}) = corr(TBL.(TBL.Properties.VariableNames{K}),TBL.(trgt),'Type','Spearman');
        end
        [temp,temp_i] = maxk(abs(cellfun(@(x)(corrResults.(season).(trgt).(x)),fieldnames(corrResults.(season).(trgt)))),num_top);
        fprintf('- %s, top %d correlated features:', trgt, num_top)
        for ii = 1:num_top
            fprintf('%s,', string(TBL.Properties.VariableNames(temp_i(ii)+1)))
        end
        fprintf('\n')
    end
    tmp = struct2table([corrResults.(season).T0, corrResults.(season).T_12, corrResults.(season).T_22, corrResults.(season).T_32, corrResults.(season).T_62, corrResults.(season).T_72]);
    AVG_corr.Season(S) = season;
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
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on','Ylim',[-1 1])
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

exportgraphics(gcf,strcat('CORRseasons_average6depths_Deadhorse.pdf'),"ContentType","image")

%% Season-based splitting for model learning

% Generate combined features and lags
clear
clc
lengths = [7*12];%[1,4,7,10,14,21,28,35,42,49];
resolution = 30;%[1,7,14,30];
res = strcat(num2str(resolution(1)),'D');
seasons = {'Spring','Summer','Autumn','Winter'};
trgts = {'T0','T_12','T_22','T_32','T_62','T_72'};
for K = 1:length(lengths)
    LEN = lengths(K);
    fprintf('**************LEN=%d**********\n',LEN)
    for S = 1:length(seasons)
        season = seasons{S};
        for J = 1:length(trgts)
            trgt = trgts{J};
            path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',"Air2UG_alltimeseries_Deadhorse_MERRA2data.csv");
            opts = detectImportOptions(path1);
            tbl = readtimetable(path1,opts);
            dt = days(1);
            tbl = retime(tbl,'regular','fillwithmissing','TimeStep',dt);
            %RESAMPLE
            tbl = movingMeanStride(tbl, str2num(cell2mat(regexp(res,'\d*','match'))));
            feature_names = ["T0","T_02","T_07","T_12","T_22","T_32","T_42","T_62","T_67","T_72","Tair"];
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
            writetimetable(tbl2,fullfile(strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\','Air2UG_alltimeseries_Deadhorse_MERRA2data_',trgt,RES,'_SeqLen', num2str(LEN), '_',season,'.csv')))
            fprintf('\n')
            clear opts dt path1
        end
    end
end







%% functions

tbl = movingMeanStride(tbl, res);

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