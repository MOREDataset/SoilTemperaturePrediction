%% Entropy/conditional entropy plots - Toolik
clear
loc = 'Toolik';
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\PreAnalysis',strcat('Entropies_ConditionalEntropies_',loc));
opts = detectImportOptions(path1);
opts.VariableNamingRule = 'preserve';
tbl = readtable(path1,opts);
tbl = movevars(tbl, "Original Entropy", "Before", "Conditional Entropy");
tbl.("Soil Temp") = categorical(tbl.("Soil Temp"));
tbl.("Feature") = categorical(tbl.("Feature"));

% compute rate of reduction for comparison reasons
tbl.("Entropy Reduction") = (tbl.("Original Entropy")-tbl.("Conditional Entropy"))./tbl.("Original Entropy")*100;

%remove SLP
% tbl(tbl.("Feature")=='SLP',:) = [];

% Assuming 'Soil Temp' and 'Feature' are categorical and 'Entropy Reduction' is numeric

% Find the unique soil temperatures and features
soilTemps = unique(tbl.('Soil Temp'));
features = categorical({'SNODP', 'SWGDN', 'T2M', 'SWLAND', 'TLML'});%, 'SLP'
features = reordercats(features, {'SNODP', 'SWGDN', 'T2M', 'SWLAND', 'TLML'});%, {'SNODP', 'SWGDN', 'T2M', 'SWLAND', 'TLML'}'SLP'

% Initialize a matrix to hold the entropy reduction data for plotting
entropyReductionData = zeros(length(soilTemps), length(features));

% Fill the matrix with the entropy reduction data
for i = 1:length(soilTemps)
    for j = 1:length(features)
        % Extract entropy reduction for each soil temp and feature
        entropyReduction = tbl.('Entropy Reduction')( ...
            tbl.('Soil Temp') == soilTemps(i) & ...
            tbl.('Feature') == features(j));
        % If we have data for this soil temp and feature, add it to the matrix
        if ~isempty(entropyReduction)
            entropyReductionData(i, j) = entropyReduction;
        end
    end
end

% Fill the matrix with the entropy reduction data
for i = 1:length(soilTemps)
    % Extract original entropy for each soil temp 
    orignalEntropy(i) = mean(tbl.('Original Entropy')(tbl.('Soil Temp') == soilTemps(i)));
end
% Create a bar graph
size1 = 35;
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile, hold all
b = bar(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}), orignalEntropy, 'FaceColor',[56 63 68]/255,'EdgeColor','none',FaceAlpha= 0.6);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}));
ylabel('Entropy','interpreter','latex','FontSize',size1);
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on')%,'Ylim',[-1 1]
% cmap = colormap(sky(length(orignalEntropy)));
% for k = 1:length(orignalEntropy)
%     b.CData(k,:) = cmap(1*k,:);
% end
nexttile, hold all
bar(categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}),entropyReductionData, 'grouped','EdgeColor','none',FaceAlpha= 0.6);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{16}$','$T_{31}$','$T_{46}$','$T_{76}$','$T_{97}$'}));
ylabel('Entropy reduction rate (\%)','interpreter','latex','FontSize',size1);
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on','Ylim',[0 47])%
lgd = legend(features,'Location','north','Interpreter','latex','Orientation','horizontal',NumColumns=2,FontSize=size1-5);
title(lgd,'RS features')
colormap(gca, sky);

% exportgraphics(t,strcat('Enrtopy_preanalysis_',loc,'.pdf'),"ContentType","image")
% close(gcf)

%% Entropy/conditional entropy plots - Deadhorse
clear
loc = 'Deadhorse';
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\PreAnalysis',strcat('Entropies_ConditionalEntropies_',loc));
opts = detectImportOptions(path1);
opts.VariableNamingRule = 'preserve';
tbl = readtable(path1,opts);
tbl = movevars(tbl, "Original Entropy", "Before", "Conditional Entropy");
tbl.("Soil Temp") = categorical(tbl.("Soil Temp"));
tbl.("Feature") = categorical(tbl.("Feature"));

% compute rate of reduction for comparison reasons
tbl.("Entropy Reduction") = (tbl.("Original Entropy")-tbl.("Conditional Entropy"))./tbl.("Original Entropy")*100;

% Assuming 'Soil Temp' and 'Feature' are categorical and 'Entropy Reduction' is numeric

% Find the unique soil temperatures and features
soilTemps = unique(tbl.('Soil Temp'));
features = categorical({'SNODP', 'SWGDN', 'T2M', 'SWLAND', 'TLML'});%, 'SLP'
features = reordercats(features, {'SNODP', 'SWGDN', 'T2M', 'SWLAND', 'TLML'});%, 'SLP'

% Initialize a matrix to hold the entropy reduction data for plotting
entropyReductionData = zeros(length(soilTemps), length(features));

% Fill the matrix with the entropy reduction data
for i = 1:length(soilTemps)
    for j = 1:length(features)
        % Extract entropy reduction for each soil temp and feature
        entropyReduction = tbl.('Entropy Reduction')( ...
            tbl.('Soil Temp') == soilTemps(i) & ...
            tbl.('Feature') == features(j));
        % If we have data for this soil temp and feature, add it to the matrix
        if ~isempty(entropyReduction)
            entropyReductionData(i, j) = entropyReduction;
        end
    end
end

% Fill the matrix with the entropy reduction data
for i = 1:length(soilTemps)
    % Extract original entropy for each soil temp 
    orignalEntropy(i) = mean(tbl.('Original Entropy')(tbl.('Soil Temp') == soilTemps(i)));
end
% Create a bar graph
size1 = 35;
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,2,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile, hold all
b = bar(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}), orignalEntropy, 'FaceColor',[56 63 68]/255,'EdgeColor','none',FaceAlpha= 0.6);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}));
ylabel('Entropy','interpreter','latex','FontSize',size1);
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on')%,'Ylim',[-1 1]
% cmap = colormap(sky(length(orignalEntropy)));
% for k = 1:length(orignalEntropy)
%     b.CData(k,:) = cmap(1*k,:);
% end
nexttile, hold all
bar(categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}),entropyReductionData, 'grouped','EdgeColor','none',FaceAlpha= 0.6);
set(gca, 'XTickLabel', categorical({'$T_{0}$','$T_{12}$','$T_{22}$','$T_{32}$','$T_{62}$','$T_{72}$'}));
ylabel('Entropy reduction rate (\%)','interpreter','latex','FontSize',size1);
set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
set(gca,'YGrid','on','YMinorGrid','on', 'Box','on','Ylim',[0 47])%,,'Ylim',[0 1.55]
lgd = legend(features,'Location','north','Interpreter','latex','Orientation','horizontal',FontSize=size1-5, NumColumns=3);
title(lgd,'RS features')
colormap(gca, sky);
% 
exportgraphics(t,strcat('Enrtopy_preanalysis_',loc,'.pdf'),"ContentType","image")
close(gcf)


%% Boxplot season-sepcific variation -Toolik
clc
clear
loc = 'Toolik';
path1 = strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\Air2UG_alltimeseries_',loc,'_MERRA2data.csv');
opts = detectImportOptions(path1);
opts.VariableNamingRule = 'preserve';
data = readtable(path1,opts);


% Convert the 'Date' column to datetime and extract 'Year', 'Month', and 'Day'
data.Date = datetime(data.Date, 'InputFormat', 'dd-MMM-yyyy');
data.Month = month(data.Date);

% Map the 'Month' to 'Season' using a switch case
data.Season = cell(height(data), 1);
for i = 1:height(data)
    month = data.Month(i);
    switch month
        case {12, 1, 2}
            data.Season{i} = 'Winter';
        case {3, 4, 5}
            data.Season{i} = 'Spring';
        case {6, 7, 8}
            data.Season{i} = 'Summer';
        case {9, 10, 11}
            data.Season{i} = 'Fall';
    end
end
data.Season = categorical(data.Season);

% Select the soil temperature columns
data2 = [data.T0;data.T16;data.T31;data.T46;data.T76;data.T97];

soildepths = strcat(strings(height(data2),1),'T0');
soildepths(height(data)+1:2*height(data)) = strcat(strings(height(data),1),'T16');
soildepths(2*height(data)+1:3*height(data)) = strcat(strings(height(data),1),'T31');
soildepths(3*height(data)+1:4*height(data)) = strcat(strings(height(data),1),'T46');
soildepths(4*height(data)+1:5*height(data)) = strcat(strings(height(data),1),'T76');
soildepths(5*height(data)+1:6*height(data)) = strcat(strings(height(data),1),'T97');
seasons = repmat(data.Season, 6, 1);

% Create the boxplot
size1 = 35;
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile, hold all
b = boxchart(seasons,data2,'GroupByColor',soildepths,'MarkerStyle','.','JitterOutliers','on','MarkerSize',20,'BoxFaceAlpha',0.6);
set(gca,'FontSize',size1,'TickLabelInterpreter','latex',ylim=[-50 30])
set(gca,'YGrid','on','YMinorGrid','on')
box("on")
axis tight

% Add labels and title
ylabel('Soil temperatures ($^{\circ}$C)','FontSize',size1-5,'Interpreter','latex')
%title('Distribution of Soil Temperatures Across Seasons')
lgd = legend({'$0cm$','$16cm$','$31cm$','$46cm$','$76cm$','$97cm$'}, 'Location','northoutside','Interpreter','latex','Orientation','horizontal',FontSize=size1-5);
title(lgd,'Soil deph')

% Create line
annotation(gcf,'line',[0.532949579831933 0.532949579831933],...
    [0.870327254305986 0.0699088145896751],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle',':');

% Create line
annotation(gcf,'line',[0.300357142857142 0.300357142857142],...
    [0.871340425531921 0.0709219858156098],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle',':');

% Create line
annotation(gcf,'line',[0.766668067226886 0.766668067226886],...
    [0.871340425531929 0.0709219858156176],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle',':');

% 
% exportgraphics(t,strcat('SoilTempDistributionSeasons_preanalysis_',loc,'.pdf'),"ContentType","image")
% close(gcf)

%% Boxplot season-sepcific variation -Deadhorse
clc
clear
loc = 'Deadhorse';
path1 = strcat('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\Air2UG_alltimeseries_',loc,'_MERRA2data.csv');
opts = detectImportOptions(path1);
opts.VariableNamingRule = 'preserve';
data = readtable(path1,opts);


% Convert the 'Date' column to datetime and extract 'Year', 'Month', and 'Day'
data.Date = datetime(data.Date, 'InputFormat', 'dd-MMM-yyyy');
data.Month = month(data.Date);

% Map the 'Month' to 'Season' using a switch case
data.Season = cell(height(data), 1);
for i = 1:height(data)
    month = data.Month(i);
    switch month
        case {12, 1, 2}
            data.Season{i} = 'Winter';
        case {3, 4, 5}
            data.Season{i} = 'Spring';
        case {6, 7, 8}
            data.Season{i} = 'Summer';
        case {9, 10, 11}
            data.Season{i} = 'Fall';
    end
end
data.Season = categorical(data.Season);

% Select the soil temperature columns
data2 = [data.T0;data.T_12;data.T_22;data.T_32;data.T_62;data.T_72];

soildepths = strcat(strings(height(data2),1),'T0');
soildepths(height(data)+1:2*height(data)) = strcat(strings(height(data),1),'T12');
soildepths(2*height(data)+1:3*height(data)) = strcat(strings(height(data),1),'T22');
soildepths(3*height(data)+1:4*height(data)) = strcat(strings(height(data),1),'T32');
soildepths(4*height(data)+1:5*height(data)) = strcat(strings(height(data),1),'T62');
soildepths(5*height(data)+1:6*height(data)) = strcat(strings(height(data),1),'T72');
seasons = repmat(data.Season, 6, 1);

% Create the boxplot
size1 = 35;
figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
nexttile, hold all
b = boxchart(seasons,data2,'GroupByColor',soildepths,'MarkerStyle','.','JitterOutliers','on','MarkerSize',20,'BoxFaceAlpha',0.6);
set(gca,'FontSize',size1,'TickLabelInterpreter','latex',ylim=[-40 20])
set(gca,'YGrid','on','YMinorGrid','on')
box("on")

% Add labels and title
ylabel('Soil temperatures ($^{\circ}$C)','FontSize',size1-5,'Interpreter','latex')
%title('Distribution of Soil Temperatures Across Seasons')
lgd = legend({'$0cm$','$12cm$','$22cm$','$32cm$','$62cm$','$72cm$'}, 'Location','northoutside','Interpreter','latex','Orientation','horizontal',FontSize=size1-5);
title(lgd,'Soil depth')


% Create line
annotation(gcf,'line',[0.532949579831933 0.532949579831933],...
    [0.870327254305986 0.0699088145896751],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle',':');

% Create line
annotation(gcf,'line',[0.300357142857142 0.300357142857142],...
    [0.871340425531921 0.0709219858156098],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle',':');

% Create line
annotation(gcf,'line',[0.766668067226886 0.766668067226886],...
    [0.871340425531929 0.0709219858156176],...
    'Color',[0.501960784313725 0.501960784313725 0.501960784313725],...
    'LineWidth',2,...
    'LineStyle',':');

% exportgraphics(t,strcat('SoilTempDistributionSeasons_preanalysis_',loc,'.pdf'),"ContentType","image")
% close(gcf)
