%% Plot observations vs predictions - Toolik lake
clear
clc
seasons = {'Winter','Spring','Summer','Autumn'};
targets = {'T0','T16','T31','T46','T76','T97'};
Lens = 0;%[0	3	6	9	13	20	27	34 41 48];
loc = 'Toolik';

size1 = 35;
mrks = {'v','square', '^', 'o'};
cmap = get(0, 'defaultaxescolororder');
cmap(2,:) = [0.4660 0.6740 0.1880];
cmap(4,:) = [0.8500 0.3250 0.0980];

%Baselines
dir_name1 = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\Baselines\';
filenames1 = {dir(fullfile(dir_name1, '*.csv')).name}.';
filenames1 = filenames1(contains(filenames1,loc));
filenames1 = filenames1(contains(filenames1(contains(filenames1,loc)),'Baseline1'));
filenames2 = {dir(fullfile(dir_name1, '*.csv')).name}.';
filenames2 = filenames2(contains(filenames2,loc));
filenames2 = filenames2(contains(filenames2,'Baseline2'));

S1 = timerange(datetime(2015,1,1),datetime(2015,12,31),'closed');
S2 = timerange(datetime(2014,1,1),datetime(2014,12,31),'closed');
%

dir_name = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\InputLength\';
files = dir(fullfile(dir_name, '*.csv'));
filenames = {files.name}.';
filenames = filenames(contains(filenames,loc));

for k = [1 3 5]
    if k==3
        S = S2;
    else
        S = S1; 
    end
    figure('units','normalized','outerposition',[0 0 1 1])
    ax1 = axes('Parent',gcf,'Units','Normalize','Position',[0.05 0.0847457627118644 0.939583333333333 0.897308075772681]);
    % t = tiledlayout('flow','TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
    Target_files = filenames(contains(filenames, targets{k})');
    Baselines_files1 = filenames1(contains(filenames1, targets{k})');
    Baselines_files2 = filenames2(contains(filenames2, targets{k})');
    ax1 = nexttile; hold on
    for j = 1:length(seasons)
        Season_files = Target_files(contains(Target_files, seasons{j})');
        Season_files = Season_files(contains(Season_files,'_1_'));

        Baselines_season1 = Baselines_files1(contains(Baselines_files1, seasons{j})');
        Baselines_season2 = Baselines_files2(contains(Baselines_files2, seasons{j})');

        path1 = strcat(dir_name,Season_files{:});
        opts = detectImportOptions(path1);
        opts.VariableNames = ["Date", "actual", "preds"];
        opts.VariableTypes = ["datetime", "double", "double"];
        opts = setvaropts(opts, "Date", "InputFormat", "yyyy.0-MM.0-dd.0");
        eval(sprintf('tbl_%s = readtimetable(path1,opts);',seasons{j}))
        eval(sprintf('tbl_%s.Season = j.*ones(height(tbl_%s),1);',seasons{j},seasons{j}))

        %
        path1_b1 = strcat(dir_name1,Baselines_season1{:});
        opts = detectImportOptions(path1_b1);
        opts.VariableNames = ["Date", "actual", "preds"];
        opts.VariableTypes = ["datetime", "double", "double"];
        opts = setvaropts(opts, "Date", "InputFormat", "yyyy.0-MM.0-dd.0");
        eval(sprintf('Baseline1_%s = readtimetable(path1_b1,opts);',seasons{j}))
        eval(sprintf('Baseline1_%s.Season = j.*ones(height(Baseline1_%s),1);',seasons{j},seasons{j}))

        %
        path1_b2 = strcat(dir_name1,Baselines_season2{:});
        opts = detectImportOptions(path1_b2);
        opts.VariableNames = ["Date", "actual", "preds"];
        opts.VariableTypes = ["datetime", "double", "double"];
        opts = setvaropts(opts, "Date", "InputFormat", "yyyy.0-MM.0-dd.0");
        eval(sprintf('Baseline2_%s = readtimetable(path1_b2,opts);',seasons{j}))
        eval(sprintf('Baseline2_%s.Season = j.*ones(height(Baseline2_%s),1);',seasons{j},seasons{j}))
    end
    tbl_allSeasons = [tbl_Winter;tbl_Spring;tbl_Summer;tbl_Autumn];
    dt = days(1);
    tbl_allSeasons = retime(tbl_allSeasons(S,:),'regular','fillwithmissing','TimeStep',dt);
    L1 = plot(tbl_allSeasons.Date,tbl_allSeasons.actual,'DisplayName','Observations','LineStyle','-','Color','k','LineWidth', 5, 'Parent',ax1);
    hold on
    
    %
    Baseline1_allSeasons = [Baseline1_Winter;Baseline1_Spring;Baseline1_Summer;Baseline1_Autumn];
    Baseline1_allSeasons = retime(Baseline1_allSeasons,'regular','fillwithmissing','TimeStep',dt);
    Baseline1_allSeasons = Baseline1_allSeasons(S,:);

    Baseline2_allSeasons = [Baseline2_Winter;Baseline2_Spring;Baseline2_Summer;Baseline2_Autumn];
    Baseline2_allSeasons = retime(Baseline2_allSeasons,'regular','fillwithmissing','TimeStep',dt);
    Baseline2_allSeasons = Baseline2_allSeasons(S,:);
    L2 = plot(Baseline1_allSeasons.Date,Baseline1_allSeasons.preds,'DisplayName','Baseline1','LineStyle',':','Marker','+','MarkerSize',10,'Color',[0.650980392156863 0.650980392156863 0.650980392156863],'LineWidth',1, 'Parent',ax1);
    L3 = plot(Baseline2_allSeasons.Date,Baseline2_allSeasons.preds,'DisplayName','Baseline2','LineStyle','--','Marker','none','MarkerSize',10,'Color',[0.650980392156863 0.650980392156863 0.650980392156863],'LineWidth',2, 'Parent',ax1);    

    for j = 1:length(seasons)
        indx = find(tbl_allSeasons.Season==j);
        tbl_tmp = retime(tbl_allSeasons(indx,1:2),'regular','fillwithmissing','TimeStep',dt);
        eval(sprintf('AM%d = plot(tbl_tmp.Date(S),tbl_tmp.preds(S),''DisplayName'',strcat(''%s'','' preds.''),''LineStyle'',''-'',''Color'',cmap(j,:),''Marker'', mrks{j},''MarkerFaceColor'',''w'',''MarkerEdgeColor'',cmap(j,:),''MarkerSize'',4,''LineWidth'',3, ''Parent'', ax1);',j,seasons{j}))
    end  
    
    set(gca,'FontSize',size1,'TickLabelInterpreter','latex');
    if k == 1
        ylabel(strcat({'Soil temperatures at '},strrep(strrep(targets{k},'_',''),'T',''),'cm',{' ($^{\circ}$C)'}),'FontSize',size1,'Interpreter','latex');
        legend1 = legend(gca,'show');%legend1 = legend(ax1,[L1 L2 L3], {'Observations','Baseline4','Baseline5',''});
        set(legend1,'Orientation','vertical','Location','northeastoutside','Interpreter','latex','FontSize',size1);
   else
        ylabel(strcat({'Soil temperatures at -'},strrep(strrep(targets{k},'_',''),'T',''),'cm',{' ($^{\circ}$C)'}),'FontSize',size1,'Interpreter','latex');
    end        

    eval(sprintf('set(gca,''Ylim'',[-40 20],''Box'',''on'',''XTickLabelRotation'',0,''FontSize'',size1-5,''TickLabelInterpreter'',''latex'',''LineWidth'',0.5,''YMinorGrid'',''on'',''YMinorTick'',''on'',''TickDir'',''in'',''TickLength'',[0.005 0.005]);'));

    %
    % ax2 = copyobj( ax1, gcf);
    % set(ax2,'Units','Normalize','Position',[0.05 0.0847457627118644 0.939583333333333 0.897308075772681],'visible','off')
    % delete( get(ax2, 'Children'))
    % if k == 3 || k == 5
    %     lh2 = legend(ax2, [AM1 AM2 AM3 AM4], strrep(seasons,'Autumn','Fall'),'Position',[0.800147847765798 0.141337367538986 0.124572285502873 0.293313090441315]);
    % else
    %     lh2 = legend(ax2, [AM1 AM2 AM3 AM4], strrep(seasons,'Autumn','Fall'),'Position',[0.801723478017898 0.575986436707381 0.124572285502873 0.293313090441315]);
    % end
    % title(lh2,'Predictions');
    % ylabel(ax1, strcat(targets{k},' ($^{\circ}C$)'),'interpreter','latex','FontSize',size1);
    % xlabel(ax1, 'Time','interpreter','latex','FontSize',size1);

    exportgraphics(gcf,strcat('PlotPredictions_Toolik_',targets{k},'.pdf'))
    % close(gcf)
end

%% Input length graph - Toolik lake - performance computation using results
clc
clear
seasons = {'Winter','Spring','Summer','Autumn'};
targets = {'T0','T16','T31','T46','T76','T97'};
Lens = [0 6 13 20 27 34 41 48 55 62 69 76 83];%[0	3	6	9	13	20	27	34 41 48];
loc = 'Toolik';

size1 = 35;
mrks = {'v','square', '^', 'o'};
cmap = get(0, 'defaultaxescolororder');
cmap(2,:) = [0.4660 0.6740 0.1880];
cmap(4,:) = [0.8500 0.3250 0.0980];

dir_name = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\InputLength\';
files = dir(fullfile(dir_name, '*.csv'));
filenames = {files.name}.';
filenames = filenames(contains(filenames,loc));

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

for k = 1:length(targets)
    Target_files = filenames(contains(filenames, targets{k})');
    nexttile, hold on
    for j = 1:length(seasons)
        Season_files = Target_files(contains(Target_files, seasons{j})');
        tmps = regexp(Season_files,'_','split');
        [LEN_value,LEN_idx] = sort(cellfun(@str2num, cellfun(@(x) x{6}, tmps, 'UniformOutput', false)));
        l = 1;
        remove_these = [];
        for i = 1:length(LEN_idx)
            if ~ismember(LEN_value(i)-1,Lens)
                remove_these(l) = i;
                l = l+1;
            end
        end
        LEN_value(remove_these) = [];
        LEN_idx(remove_these) = [];
       
        for i = 1:length(LEN_idx)
            path1 = strcat(dir_name,Season_files{LEN_idx(i)});
            opts = detectImportOptions(path1);
            opts.VariableNames = ["Date", "actual", "preds"];
            opts.VariableTypes = ["datetime", "double", "double"];
            opts = setvaropts(opts, "Date", "InputFormat", "yyyy.0-MM.0-dd.0");
            tbl = readtimetable(path1,opts);
            delta = tbl.actual-tbl.preds;
            RMSE(i, j) = sqrt(mean((delta).^2)); %row:LEN, column:season
            R2 = 1-sum((delta).^2)./sum((tbl.actual-mean(tbl.actual)).^2);
            % percentile25(i, j) = mean(delta(sign(delta)==-1));
            % percentile75(i,j) = mean(delta(sign(delta)==1));
        end
        
        % findout the optimal look-back period index
        tmp = find_optimal(RMSE, j, 0.1, 0.03);
        
        %plot
        plot(categorical(Lens),RMSE(:,j)','DisplayName',strrep(seasons{j},'Autumn','Fall'),'LineStyle','-','Color',cmap(j,:),'Marker', mrks{j},'MarkerFaceColor','w','MarkerEdgeColor',cmap(j,:),'MarkerSize',7,'LineWidth', 2);
        plot([categorical(Lens(tmp)), categorical(Lens(tmp))], [0, RMSE(tmp,j)], ':','LineWidth', 1.2,'Color',cmap(j,:),'HandleVisibility', 'Off')
        scatter(categorical(Lens(tmp)),RMSE(tmp,j),70-j*10,"filled",'Marker', mrks{j},'MarkerEdgeColor',cmap(j,:),'MarkerFaceColor',cmap(j,:),'HandleVisibility', 'Off')
        % errorbar(categorical(Lens),RMSE(:,j)',percentile25(:,j)', percentile75(:,j)','LineStyle','none','Color',cmap(j,:),'HandleVisibility', 'Off')
    end
    set(gca, 'XDir','reverse')
    set(gca,'FontSize',size1,'TickLabelInterpreter','latex');%,'Ylim',[0 7]); just in the case of Household 3
    if k == 2 
        legend1 = legend(gca,'show');
        set(legend1,'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',size1);
    end
    eval(sprintf('set(gca,''Box'',''on'',''XTickLabelRotation'',0,''FontSize'',size1-5,''TickLabelInterpreter'',''latex'',''LineWidth'',0.5,''YMinorGrid'',''on'',''YMinorTick'',''on'',''TickDir'',''in'',''TickLength'',[0.005 0.005]);'));
    eval(sprintf('title(gca,''%s'',''Interpreter'',''latex'',''FontSize'',size1);', targets{k}))
end
ylabel(t, 'RMSE ($^{\circ}C$)','interpreter','latex','FontSize',size1);
xlabel(t, 'Look-back period (days)','interpreter','latex','FontSize',size1);

exportgraphics(t,strcat('InputLength_',loc,'.pdf'))
close(gcf)

%% Training set size - Toolik lake
clear
seasons = {'Winter','Spring','Summer','Autumn'};
targets = {'T0','T16','T31', 'T46', 'T76', 'T97'};
loc = 'Toolik';

size1 = 35;
mrks = {'v','square', '^', 'o'};
cmap = get(0, 'defaultaxescolororder');
cmap(2,:) = [0.4660 0.6740 0.1880];
cmap(4,:) = [0.8500 0.3250 0.0980];

dir_name = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\TrainSize\';
files = dir(fullfile(dir_name, '*.csv'));
filenames = {files.name}.';
filenames = filenames(contains(filenames,loc));

cmp = [
    % 0, 0, 0;        % Black
    % 144, 0, 99;     % Dark purple
    % 227, 0, 119;    % Purple-red
    252, 129, 89;   % Red-orange
    253, 199, 79;   % Yellow-orange
    251, 252, 191   % Bright yellow
    173, 216, 230;    % Light blue (sky color, for a yellow-blue transition)
    0, 128, 255       % Bright blue
    ] / 255; % Normalize the RGB values to be between 0 and 1
cmap = interp1(linspace(1,25,size(cmp,1)),flip(cmp),linspace(1,25,25));


figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('flow','TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

for j = 1:length(seasons)
    Season_files = filenames(contains(filenames, seasons{j})');
    eval(sprintf('ax%d = nexttile; hold on',j))
    for k = 1:length(targets)
        Target_files = Season_files(contains(Season_files, targets{k})');
        tmps = str2double(string(regexp(string(regexp(Target_files,'\d*years','match')),'\d*','match')));          
        [LEN_value,LEN_idx] = sort(tmps);      
        for i = 1:length(LEN_idx)
            path1 = strcat(dir_name,Target_files{LEN_idx(i)});
            opts = detectImportOptions(path1);
            opts.VariableNames = ["Date", "actual", "preds"];
            opts.VariableTypes = ["datetime", "double", "double"];
            opts = setvaropts(opts, "Date", "InputFormat", "yyyy.0-MM.0-dd.0");
            tbl = readtimetable(path1,opts);
            delta = tbl.actual-tbl.preds;
            RMSE(i, k, j) = sqrt(mean((delta).^2)); %row:LEN, column:season
            % percentile25(i, j) = mean(delta(sign(delta)==-1));
            % percentile75(i,j) = mean(delta(sign(delta)==1));
        end
        if LEN_idx < 14
            RMSE(14, :, j) = NaN(1,length(targets));
        end
    end
    surf(categorical(targets), categorical(1:14), RMSE(:,:,j), 'FaceColor','interp')
    view(2)
    colormap(cmap);
    hcb = colorbar;
    hcb.Label.String = 'RMSE ($^{\circ}C$)';
    hcb.Label.Interpreter = 'latex';

    % xticklabels(categorical({'T0','','T16','','T31','','T46','','T76','','T97'}));
    yticklabels(categorical(1:14));
    title(strrep(seasons{j},'Autumn','Fall'),'Interpreter','latex')
    set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
    eval(sprintf('set(gca,''Box'',''on'',''XTickLabelRotation'',0,''FontSize'',size1-10,''TickLabelInterpreter'',''latex'',''LineWidth'',0.5,''TickDir'',''in'',''TickLength'',[0.005 0.005]);'));
end
for j = 1:length(seasons)
    eval(sprintf('clim(ax%d, [0, round(max(RMSE(:,:,:),[],"all"))]);',j))
end

xlabel(t, 'Soil temperature','interpreter','latex','FontSize',size1);
ylabel(t, 'Training size (years)','interpreter','latex','FontSize',size1);

exportgraphics(t,strcat('TrainSize_',loc,'.pdf'),"ContentType","image")
close(gcf)

%% Plot observations vs predictions - Deadhorse
clear
clc
seasons = {'Winter','Spring','Summer','Autumn'};
targets = {'T0','T_12','T_22', 'T_32', 'T_62', 'T_72'};
Lens = 0;%[0	3	6	9	13	20	27	34 41 48];
loc = 'Deadhorse';

size1 = 35;
mrks = {'v','square', '^', 'o'};
cmap = get(0, 'defaultaxescolororder');
cmap(2,:) = [0.4660 0.6740 0.1880];
cmap(4,:) = [0.8500 0.3250 0.0980];

%Baselines
dir_name1 = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\Baselines\';
filenames1 = {dir(fullfile(dir_name1, '*.csv')).name}.';
filenames1 = filenames1(contains(filenames1,loc));
filenames1 = filenames1(contains(filenames1(contains(filenames1,loc)),'Baseline1'));
filenames2 = {dir(fullfile(dir_name1, '*.csv')).name}.';
filenames2 = filenames2(contains(filenames2,loc));
filenames2 = filenames2(contains(filenames2,'Baseline2'));

S1 = timerange(datetime(1999,1,1),datetime(1999,12,31),'closed');
S2 = timerange(datetime(1998,1,1),datetime(1998,12,31),'closed');
%

dir_name = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\InputLength\';
files = dir(fullfile(dir_name, '*.csv'));
filenames = {files.name}.';
filenames = filenames(contains(filenames,loc));

for k = [1 4 6]
    if k==3
        S = S2;
    else
        S = S1; 
    end
    figure('units','normalized','outerposition',[0 0 1 1])
    ax1 = axes('Parent',gcf,'Units','Normalize','Position',[0.05 0.0847457627118644 0.939583333333333 0.897308075772681]);
    % t = tiledlayout('flow','TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);
    Target_files = filenames(contains(filenames, targets{k})');
    Baselines_files1 = filenames1(contains(filenames1, targets{k})');
    Baselines_files2 = filenames2(contains(filenames2, targets{k})');
    ax1 = nexttile; hold on
    for j = 1:length(seasons)
        Season_files = Target_files(contains(Target_files, seasons{j})');
        Season_files = Season_files(contains(Season_files,'_1_'));

        Baselines_season1 = Baselines_files1(contains(Baselines_files1, seasons{j})');
        Baselines_season2 = Baselines_files2(contains(Baselines_files2, seasons{j})');

        path1 = strcat(dir_name,Season_files{:});
        opts = detectImportOptions(path1);
        opts.VariableNames = ["Date", "actual", "preds"];
        opts.VariableTypes = ["datetime", "double", "double"];
        opts = setvaropts(opts, "Date", "InputFormat", "yyyy.0-MM.0-dd.0");
        eval(sprintf('tbl_%s = readtimetable(path1,opts);',seasons{j}))
        eval(sprintf('tbl_%s.Season = j.*ones(height(tbl_%s),1);',seasons{j},seasons{j}))

        %
        path1_b1 = strcat(dir_name1,Baselines_season1{:});
        opts = detectImportOptions(path1_b1);
        opts.VariableNames = ["Date", "actual", "preds"];
        opts.VariableTypes = ["datetime", "double", "double"];
        opts = setvaropts(opts, "Date", "InputFormat", "yyyy-MM-dd");
        eval(sprintf('Baseline1_%s = readtimetable(path1_b1,opts);',seasons{j}))
        eval(sprintf('Baseline1_%s.Season = j.*ones(height(Baseline1_%s),1);',seasons{j},seasons{j}))

        %
        path1_b2 = strcat(dir_name1,Baselines_season2{:});
        opts = detectImportOptions(path1_b2);
        opts.VariableNames = ["Date", "actual", "preds"];
        opts.VariableTypes = ["datetime", "double", "double"];
        opts = setvaropts(opts, "Date", "InputFormat", "yyyy-MM-dd");
        eval(sprintf('Baseline2_%s = readtimetable(path1_b2,opts);',seasons{j}))
        eval(sprintf('Baseline2_%s.Season = j.*ones(height(Baseline2_%s),1);',seasons{j},seasons{j}))
    end
    tbl_allSeasons = [tbl_Winter;tbl_Spring;tbl_Summer;tbl_Autumn];
    dt = days(1);
    tbl_allSeasons = retime(tbl_allSeasons(S,:),'regular','fillwithmissing','TimeStep',dt);
    L1 = plot(tbl_allSeasons.Date,tbl_allSeasons.actual,'DisplayName','Observations','LineStyle','-','Color','k','LineWidth', 5, 'Parent',ax1);
    hold on
    
    %
    Baseline1_allSeasons = [Baseline1_Winter;Baseline1_Spring;Baseline1_Summer;Baseline1_Autumn];
    Baseline1_allSeasons = retime(Baseline1_allSeasons,'regular','fillwithmissing','TimeStep',dt);
    Baseline1_allSeasons = Baseline1_allSeasons(S,:);

    Baseline2_allSeasons = [Baseline2_Winter;Baseline2_Spring;Baseline2_Summer;Baseline2_Autumn];
    Baseline2_allSeasons = retime(Baseline2_allSeasons,'regular','fillwithmissing','TimeStep',dt);
    Baseline2_allSeasons = Baseline2_allSeasons(S,:);
    L2 = plot(Baseline1_allSeasons.Date,Baseline1_allSeasons.preds,'DisplayName','Baseline1','LineStyle',':','Marker','+','MarkerSize',10,'Color',[0.650980392156863 0.650980392156863 0.650980392156863],'LineWidth',1, 'Parent',ax1);
    L3 = plot(Baseline2_allSeasons.Date,Baseline2_allSeasons.preds,'DisplayName','Baseline2','LineStyle','--','Marker','none','MarkerSize',10,'Color',[0.650980392156863 0.650980392156863 0.650980392156863],'LineWidth',2, 'Parent',ax1);    

    for j = 1:length(seasons)
        indx = find(tbl_allSeasons.Season==j);
        tbl_tmp = retime(tbl_allSeasons(indx,1:2),'regular','fillwithmissing','TimeStep',dt);
        eval(sprintf('AM%d = plot(tbl_tmp.Date(S),tbl_tmp.preds(S),''DisplayName'',strcat(''%s'','' preds.''),''LineStyle'',''-'',''Color'',cmap(j,:),''Marker'', mrks{j},''MarkerFaceColor'',''w'',''MarkerEdgeColor'',cmap(j,:),''MarkerSize'',4,''LineWidth'',3, ''Parent'', ax1);',j,seasons{j}))
    end  
    
    set(gca,'FontSize',size1,'TickLabelInterpreter','latex');
    if k == 1
        ylabel(strcat({'Soil temperatures at '},strrep(strrep(targets{k},'_',''),'T',''),'cm',{' ($^{\circ}$C)'}),'FontSize',size1,'Interpreter','latex');
        legend1 = legend(gca,'show');%legend1 = legend(ax1,[L1 L2 L3], {'Observations','Baseline4','Baseline5',''});
        set(legend1,'Orientation','vertical','Location','northeastoutside','Interpreter','latex','FontSize',size1);
   else
        ylabel(strcat({'Soil temperatures at -'},strrep(strrep(targets{k},'_',''),'T',''),'cm',{' ($^{\circ}$C)'}),'FontSize',size1,'Interpreter','latex');
    end        

    eval(sprintf('set(gca,''Ylim'',[-40 20],''Box'',''on'',''XTickLabelRotation'',0,''FontSize'',size1-5,''TickLabelInterpreter'',''latex'',''LineWidth'',0.5,''YMinorGrid'',''on'',''YMinorTick'',''on'',''TickDir'',''in'',''TickLength'',[0.005 0.005]);'));
    % 
    % %
    % ax2 = copyobj( ax1, gcf);
    % set(ax2,'Units','Normalize','Position',[0.05 0.0847457627118644 0.939583333333333 0.897308075772681],'visible','off')
    % delete( get(ax2, 'Children'))
    % % if k ~= 1
    % lh2 = legend(ax2, [AM1 AM2 AM3 AM4], strrep(seasons,'Autumn','Fall'),'Position',[0.800147847765798 0.141337367538986 0.124572285502873 0.293313090441315]);
    % % else
    % % lh2 = legend(ax2, [AM1 AM2 AM3 AM4], strrep(seasons,'Autumn','Fall'),'Position',[0.801723478017898 0.575986436707381 0.124572285502873 0.293313090441315]);
    % % end
    % title(lh2,'Predictions');
    % ylabel(ax1, strcat(strrep(targets{k},'_',''),' ($^{\circ}C$)'),'interpreter','latex','FontSize',size1);
    % xlabel(ax1, 'Time','interpreter','latex','FontSize',size1);

    exportgraphics(gcf,strcat('PlotPredictions_Deadhorse_',strrep(targets{k},'_',''),'.pdf'))
    close(gcf)
end


%% Input length graph - Deadhorse - performance computation using results
clear
seasons = {'Winter','Spring','Summer','Autumn'};
targets = {'T0','T_12','T_22', 'T_32', 'T_62', 'T_72'};
Lens = [0 6 13 20 27 34 41 48 55 62 69 76 83];%[0	3	6	9	13	20	27	34 41 48];
loc = 'Deadhorse';

size1 = 35;
mrks = {'v','square', '^', 'o'};
cmap = get(0, 'defaultaxescolororder');
cmap(2,:) = [0.4660 0.6740 0.1880];
cmap(4,:) = [0.8500 0.3250 0.0980];

dir_name = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\InputLength\';
files = dir(fullfile(dir_name, '*.csv'));
filenames = {files.name}.';
filenames = filenames(contains(filenames,loc));

figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

for k = 1:length(targets)
    Target_files = filenames(contains(filenames, targets{k})');
    nexttile, hold on
    for j = 1:length(seasons)
        Season_files = Target_files(contains(Target_files, seasons{j})');
        tmps = regexp(Season_files,'_','split');
        if strcmp(targets{k},'T0')
            [LEN_value,LEN_idx] = sort(cellfun(@str2num, cellfun(@(x) x{6}, tmps, 'UniformOutput', false)));
        else
            [LEN_value,LEN_idx] = sort(cellfun(@str2num, cellfun(@(x) x{7}, tmps, 'UniformOutput', false)));
        end            
        l = 1;
        remove_these = [];
        for i = 1:length(LEN_idx)
            if ~ismember(LEN_value(i)-1,Lens)
                remove_these(l) = i;
                l = l+1;
            end
        end
        LEN_value(remove_these) = [];
        LEN_idx(remove_these) = [];
       
        for i = 1:length(LEN_idx)
            path1 = strcat(dir_name,Season_files{LEN_idx(i)});
            opts = detectImportOptions(path1);
            opts.VariableNames = ["Date", "actual", "preds"];
            opts.VariableTypes = ["datetime", "double", "double"];
            opts = setvaropts(opts, "Date", "InputFormat", "yyyy.0-MM.0-dd.0");
            tbl = readtimetable(path1,opts);
            delta = tbl.actual-tbl.preds;
            RMSE(i, j) = sqrt(mean((delta).^2)); %row:LEN, column:season
            % percentile25(i, j) = mean(delta(sign(delta)==-1));
            % percentile75(i,j) = mean(delta(sign(delta)==1));
        end

        % findout the optimal look-back period index
        tmp = find_optimal(RMSE, j, 0.1, 0.03);

        %plot
        plot(categorical(Lens),RMSE(:,j)','DisplayName',strrep(seasons{j},'Autumn','Fall'),'LineStyle','-','Color',cmap(j,:),'Marker', mrks{j},'MarkerFaceColor','w','MarkerEdgeColor',cmap(j,:),'MarkerSize',7,'LineWidth', 2);
        plot([categorical(Lens(tmp)), categorical(Lens(tmp))], [0, RMSE(tmp,j)], ':','LineWidth', 1.2,'Color',cmap(j,:),'HandleVisibility', 'Off')
        scatter(categorical(Lens(tmp)),RMSE(tmp,j),70-j*10,"filled",'Marker', mrks{j},'MarkerEdgeColor',cmap(j,:),'MarkerFaceColor',cmap(j,:),'HandleVisibility', 'Off')
        % errorbar(categorical(Lens),RMSE(:,j)',percentile25(:,j)', percentile75(:,j)','LineStyle','none','Color',cmap(j,:),'HandleVisibility', 'Off')
    end
    set(gca, 'XDir','reverse')
    set(gca,'FontSize',size1,'TickLabelInterpreter','latex');%,'Ylim',[0 7]); just in the case of Household 3
    if k == 2 
        legend1 = legend(gca,'show');
        set(legend1,'Orientation','horizontal','Location','northoutside','Interpreter','latex','FontSize',size1);
    end
    eval(sprintf('set(gca,''Box'',''on'',''XTickLabelRotation'',0,''FontSize'',size1-5,''TickLabelInterpreter'',''latex'',''LineWidth'',0.5,''YMinorGrid'',''on'',''YMinorTick'',''on'',''TickDir'',''in'',''TickLength'',[0.005 0.005]);'));
    eval(sprintf('title(gca,''%s'',''Interpreter'',''latex'',''FontSize'',size1);', strrep(targets{k},'_','')))
end
ylabel(t, 'RMSE ($^{\circ}C$)','interpreter','latex','FontSize',size1);
xlabel(t, 'Look-back period (days)','interpreter','latex','FontSize',size1);

exportgraphics(t,strcat('InputLength_',loc,'.pdf'))
close(gcf)


%% Training set size - Deadhorse
clear
seasons = {'Winter','Spring','Summer','Autumn'};
targets = {'T0','T_12','T_22', 'T_32', 'T_62', 'T_72'};

loc = 'Deadhorse';

size1 = 35;
mrks = {'v','square', '^', 'o'};
cmap = get(0, 'defaultaxescolororder');
cmap(2,:) = [0.4660 0.6740 0.1880];
cmap(4,:) = [0.8500 0.3250 0.0980];

dir_name = 'C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures\results4\TrainSize\';
files = dir(fullfile(dir_name, '*.csv'));
filenames = {files.name}.';
filenames = filenames(contains(filenames,loc));

cmp = [
    % 0, 0, 0;        % Black
    % 144, 0, 99;     % Dark purple
    % 227, 0, 119;    % Purple-red
    252, 129, 89;   % Red-orange
    253, 199, 79;   % Yellow-orange
    251, 252, 191   % Bright yellow
    173, 216, 230;    % Light blue (sky color, for a yellow-blue transition)
    0, 128, 255       % Bright blue
    ] / 255; % Normalize the RGB values to be between 0 and 1
cmap = interp1(linspace(1,25,size(cmp,1)),flip(cmp),linspace(1,25,25));


figure('units','normalized','outerposition',[0 0 1 1])
t = tiledlayout('flow','TileSpacing','tight','Padding','tight','units','normalized','outerposition',[0 0 1 1]);

for j = 1:length(seasons)
    Season_files = filenames(contains(filenames, seasons{j})');
    eval(sprintf('ax%d = nexttile; hold on',j))
    for k = 1:length(targets)
        Target_files = Season_files(contains(Season_files, targets{k})');
        tmps = str2double(string(regexp(string(regexp(Target_files,'\d*years','match')),'\d*','match')));          
        [LEN_value,LEN_idx] = sort(tmps);      
        for i = 1:length(LEN_idx)
            path1 = strcat(dir_name,Target_files{LEN_idx(i)});
            opts = detectImportOptions(path1);
            opts.VariableNames = ["Date", "actual", "preds"];
            opts.VariableTypes = ["datetime", "double", "double"];
            opts = setvaropts(opts, "Date", "InputFormat", "yyyy.0-MM.0-dd.0");
            tbl = readtimetable(path1,opts);
            delta = tbl.actual-tbl.preds;
            RMSE(i, k, j) = sqrt(mean((delta).^2)); %row:LEN, column:season
            % percentile25(i, j) = mean(delta(sign(delta)==-1));
            % percentile75(i,j) = mean(delta(sign(delta)==1));
        end
        if LEN_idx < 12
            RMSE(12, :, j) = NaN(1,length(targets));
        end
    end
    surf(categorical(1:length(targets)), categorical(1:12), RMSE(:,:,j), 'FaceColor','interp')
    view(2)
    colormap(cmap);
    hcb = colorbar;
    hcb.Label.String = 'RMSE ($^{\circ}C$)';
    hcb.Label.Interpreter = 'latex';

    xticklabels(categorical({'T0','T12','T22','T32','T62','T72'}));
    yticklabels(categorical(1:14));
    title(strrep(seasons{j},'Autumn','Fall'),'Interpreter','latex')
    set(gca,'FontSize',size1,'TickLabelInterpreter','latex')
    eval(sprintf('set(gca,''Box'',''on'',''XTickLabelRotation'',0,''FontSize'',size1-10,''TickLabelInterpreter'',''latex'',''LineWidth'',0.5,''TickDir'',''in'',''TickLength'',[0.005 0.005]);'));
end
for j = 1:length(seasons)
    eval(sprintf('clim(ax%d, [0, round(max(RMSE(:,:,:),[],"all"))]);',j))
end

xlabel(t, 'Soil temperature','interpreter','latex','FontSize',size1);
ylabel(t, 'Training size (years)','interpreter','latex','FontSize',size1);

exportgraphics(t,strcat('TrainSize_',loc,'.pdf'),"ContentType","image")
close(gcf)

%% Functions
function optimal_index = find_optimal(RMSE, j,rmse_threshold, rmse_threshold2)
    % Initialize the index of the optimal look-back period
    optimal_index = 1; % Start with the shortest look-back period
    
    % Assume RMSE is a matrix where each column corresponds to a different model
    % and each row corresponds to a different look-back period.
    % Initialize a variable to store the best RMSE found
    best_rmse = RMSE(optimal_index, j);
    best_index = optimal_index;
    
    % Loop through the RMSE values to find the optimal look-back period
    for O = 2:length(RMSE(:, j))
        % Calculate the difference between the current RMSE and the shortest period RMSE
        rmse_diff = RMSE(1, j) - RMSE(O, j);
        % Check if the RMSE improvement is significant enough
        if rmse_diff > rmse_threshold
            % Update the best RMSE and index if the improvement is larger than the threshold
            if RMSE(O, j) < best_rmse
                best_rmse = RMSE(O, j);
                best_index = O;
            end
        end
    end
    
    % Now, find the shortest period that is within an acceptable range of the best RMSE
    for O = 1:best_index
        if abs(RMSE(O, j) - best_rmse) <= rmse_threshold2*O
            if RMSE(1, j) - RMSE(O, j) > rmse_threshold2*O 
                optimal_index = O;
            else
                optimal_index = 1;
            end
            break;
        end
    end
end