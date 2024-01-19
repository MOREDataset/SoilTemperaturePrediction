clear
clc
seasons = {'Winter','Spring','Summer','Fall'};

loc = 'Deadhorse';
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',strcat('Air2UG_alltimeseries_',loc,'_MERRA2data.csv'));
opts = detectImportOptions(path1);
tbl_deadhorse = readtimetable(path1,opts);

loc = 'Toolik';
path1 = fullfile('C:\Users\mohamed.ahajjam\Desktop\UND\Defense resiliency platform\Datasets\AirAndActiveLayerTemperatures',strcat('Air2UG_alltimeseries_',loc,'_MERRA2data.csv'));
opts = detectImportOptions(path1);
tbl_toolik = readtimetable(path1,opts);


%% Say-day air temperature baseline - Deadhorse
clc
loc = 'Deadhorse';
trgts = {'T0','T_12','T_22','T_32','T_62','T_72'};
tbl_deadhorse = Date2Season(tbl_deadhorse);

yearCount = 3; %past years to average from
sizeTrain = floor(height(tbl_deadhorse)*0.8);
Test = tbl_deadhorse(floor(sizeTrain)+1:end,:);
sameDays_3Yaverage = Compute_sameDaysPreds(tbl_deadhorse, Test, yearCount);

%Performance
fprintf('%s Same-day %dyears average prediction: Target: Season: RMSE,MSE,VAR,R2:\n',loc,yearCount)
for h = 1:length(seasons)
    for i = 1:length(trgts)
        eval(sprintf('y_test = (tbl_deadhorse.%s(Test.Season==seasons{h}));',trgts{i}))
        y_preds = sameDays_3Yaverage.(seasons{h})';
        RMSE = sqrt(mean(mean((y_test-y_preds).^2)));
        MAE = mean(abs(y_test-y_preds));
        MSE = mean((y_test-y_preds).^2);
        VAR = 1-(var(y_test-y_preds))/(var(y_test));
        R2 = mean(1-sum((y_test-y_preds).^2)./sum((y_test-mean(y_test)).^2));
        fprintf(':%s:%s: %.3f,%.3f\n',seasons{h},trgts{i},RMSE,MAE)
    end
    fprintf('\n')   
end

%% Say-day air temperature baseline - toolik
clc
loc = 'Toolik';
trgts = {'T0','T16','T31','T46','T76','T97'};
tbl_toolik = Date2Season(tbl_toolik);
%%
yearCount = 3; %past years to average from
sizeTrain = floor(height(tbl_toolik)*0.8);
Test = tbl_toolik(floor(sizeTrain)+1:end,:);
sameDays_3Yaverage = Compute_sameDaysPreds(tbl_toolik, Test, yearCount);

%Performance
fprintf('Same-day %dyears average prediction: Target: Season: RMSE,MSE,VAR,R2:\n',yearCount)
for h = 1:length(seasons)
    for i = 1:length(trgts)
        eval(sprintf('y_test = (tbl_toolik.%s(Test.Season==seasons{h}));',trgts{i}))
        y_preds = sameDays_3Yaverage.(seasons{h})';
        RMSE = sqrt(mean(mean((y_test-y_preds).^2)));
        MAE = mean(abs(y_test-y_preds));
        MSE = mean((y_test-y_preds).^2);
        VAR = 1-(var(y_test-y_preds))/(var(y_test));
        R2 = mean(1-sum((y_test-y_preds).^2)./sum((y_test-mean(y_test)).^2));
        fprintf(':%s:%s: %.3f,%.3f\n',seasons{h},trgts{i},RMSE,MAE)
    end
    fprintf('\n')
end

%% Same-day different Location - Deadhorse
trgts1 = {'T0','T_12','T_32','T_72'};
trgts2 = {'T0','T16','T31','T76'};

sizeTrain1 = floor(height(tbl_deadhorse)*0.8);
Test1 = tbl_deadhorse(floor(sizeTrain1)+1:end,:);

preds = Compute_DeadhorsefromToolikPreds(Test1, tbl_toolik, trgts2);


%Performance
fprintf('Same-day different location prediction: Target: Season: RMSE,MAE:\n')
for h = 1:length(seasons)
    y_preds = preds.(seasons{h});
    for i = 1:length(trgts1)
        eval(sprintf('y_test = (tbl_deadhorse.%s(Test1.Season==seasons{h}));',trgts1{i}))
        RMSE = sqrt(mean(mean((y_test-y_preds(:,i)).^2)));
        MAE = mean(abs(y_test-y_preds(:,i)));
        fprintf(':%s:%s: %.3f,%.3f\n',seasons{h},trgts1{i},RMSE,MAE)
    end
    fprintf('\n')
end



%% Same-day different Location - Toolik
trgts1 = {'T0','T_12','T_32','T_72'};
trgts2 = {'T0','T16','T31','T76'};

sizeTrain2 = floor(height(tbl_toolik)*0.8);
Test2 = tbl_toolik(floor(sizeTrain2)+1:end,:);

preds = Compute_ToolikfromDeadhorsePreds(Test1, tbl_toolik, trgts2);

%Performance
fprintf('Same-day different location prediction: Target: Season: RMSE,MAE:\n')
for h = 1:length(seasons)
    y_preds = preds.(seasons{h});
    for i = 1:length(trgts2)
        eval(sprintf('y_test = (tbl_toolik.%s(Test1.Season==seasons{h}));',trgts2{i}))
        RMSE = sqrt(mean(mean((y_test-y_preds(:,i)).^2)));
        MAE = mean(abs(y_test-y_preds(:,i)));
        fprintf(':%s:%s: %.3f,%.3f\n',seasons{h},trgts2{i},RMSE,MAE)
    end
    fprintf('\n')
end
%% Functions
function tbl = Date2Season(tbl)
    tbl.Month = month(tbl.Date);
    % Map the 'Month' to 'Season' using a switch case
    tbl.Season = cell(height(tbl), 1);
    for i = 1:height(tbl)
        mnth = tbl.Month(i);
        switch mnth
            case {12, 1, 2}
                tbl.Season{i} = 'Winter';
            case {3, 4, 5}
                tbl.Season{i} = 'Spring';
            case {6, 7, 8}
                tbl.Season{i} = 'Summer';
            case {9, 10, 11}
                tbl.Season{i} = 'Fall';
        end
    end
    tbl.Season = categorical(tbl.Season);
end

function sameDays_3Yaverage = Compute_sameDaysPreds(tbl_deadhorse, Test, yearCount)
    seasons = {'Winter','Spring','Summer','Fall'};
    for h = 1:length(seasons)
        indices = find(Test.Season==seasons{h});
        for i = 1:length(indices)
            t = Test.Properties.RowTimes(indices(i));
            for j = 1:yearCount
                F1(1,j) = tbl_deadhorse{timerange(t-years(j)+days(1),'days'),'Tair'};
            end
            sameDays_3Yaverage.(seasons{h})(i) = mean(F1);
        end
    end
end


function preds = Compute_DeadhorsefromToolikPreds(Test1, tbl_toolik, trgts2)
    seasons = {'Winter','Spring','Summer','Fall'};
    for h = 1:length(seasons)
        indices = find(Test1.Season==seasons{h});
        mm_dds = datestr(tbl_toolik.Date(tbl_toolik.Season==seasons{h}), 'mm-dd');
        for i = 1:length(indices)
            t = Test1.Properties.RowTimes(indices(i));
            indices2 = matches(string(mm_dds),datestr(t, 'mm-dd'));
            F1 = tbl_toolik{indices2,trgts2};
            preds.(seasons{h})(i,:) = mean(F1);
        end
    end
end


function preds = Compute_ToolikfromDeadhorsePreds(Test2, tbl_deadhorse, trgts1)
    seasons = {'Winter','Spring','Summer','Fall'};
    for h = 1:length(seasons)
        indices = find(Test2.Season==seasons{h});
        mm_dds = datestr(tbl_deadhorse.Date(tbl_deadhorse.Season==seasons{h}), 'mm-dd');
        for i = 1:length(indices)
            t = Test2.Properties.RowTimes(indices(i));
            indices2 = matches(string(mm_dds),datestr(t, 'mm-dd'));
            F1 = tbl_deadhorse{indices2,trgts1};
            preds.(seasons{h})(i,:) = mean(F1);
        end
    end
end