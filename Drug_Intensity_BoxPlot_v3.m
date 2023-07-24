% V2: for multi group plotting & statistic exam in sub functions

close all;

keyword = 'BR_AP4_*.csv';    % apparatus name

Write2CSV = 0;

listing = dir(keyword);
name = string(extractfield(listing, 'name'));
fprintf('Data in file:\n');
name_seq = [[1:length(name)]', name'];
disp(name_seq);
% close 
for ii = 1:length(name)
    
    %% File input
    TableName = name(ii);
    T = readtable(TableName,'NumHeaderLines',1);
    %     T1 = readmatrix(TableName,'NumHeaderLines',1);
    T = table2array(T);
    
% %     % Look into different FOV 
%     groupIndex = T(:,1);
%     Index = 3;
%     CD_Cy3 = T(groupIndex == Index,3);
%     CH_Cy3 = T(groupIndex == Index,4);
%     off_Cy3 = T(groupIndex == Index,5);
% %     CD_Cy5 = T(groupIndex == Index,6);
% %     CH_Cy5 = T(groupIndex == Index,7);
% %     off_Cy5 = T(groupIndex == Index,8);
%     CD_allBac = T(groupIndex == Index,6);
%     CH_allBac = T(groupIndex == Index,7);
%     off_allBac = T(groupIndex == Index,8);
    
    CH_mask = T(:,3);
    Cy3_mask = T(:,4);
    Cy5_mask = T(:,5);
    
%     CD_Cy5 = T(:,6);
%     CH_Cy5 = T(:,7);
%     off_Cy5 = T(:,8);
    
    
   % use neg ctrl sample 4 to calibrate the compensation coeff
   RescaleCoeff0 = 1;
   Cy3_mask = Cy3_mask*RescaleCoeff0;
    
    
   RescaleCoeff1 = 1;
   Cy5_mask = Cy5_mask*RescaleCoeff1;
    
   RescaleCoeff2 = 1;
   CH_mask = CH_mask*RescaleCoeff2;
    
%     %FA4
%     CD_percent_Cy3 = StatisticBoxPlot(CD_Cy3, CH_Cy3, CH_low_Cy3,2,1);
%     CD_percent_Cy5 = StatisticBoxPlot(CD_Cy5, CH_Cy5, CH_low_Cy5,1);
%     CD_percent_allBac = StatisticBoxPlot(CD_allBac, CH_allBac, CH_low_allBac,2,1);
    
    % normal
    Drug_low_Cy3 = 0;
    Cy3_output = StatisticBoxPlot_v2(Cy3_mask, Drug_low_Cy3,1,2);
    pause(3);
    
    Drug_low_Cy5 = 0;
    Cy5_output = StatisticBoxPlot_v2(Cy5_mask, Drug_low_Cy5,1,2);
    pause(3);
    
    Drug_low_CH = 0;
    CH_output = StatisticBoxPlot_v2(CH_mask, Drug_low_CH,1,2);
    
    %% Save to original csv file
    if Write2CSV == 1
    Tupdate = array2table(T);
    Tupdate(:,6) = array2table(Cy3_output);
    Tupdate(:,7) = array2table(Cy5_output);
    Tupdate(:,8) = array2table(CH_output);
    writetable(Tupdate, TableName);
    end
    
    %% Plotting default settings
    
    set(0, 'DefaultTextFontSize', 18);
    set(0, 'DefaultLineLineWidth', 2);
    set(0, 'DefaultAxesLIneWidth', 2);
    set(0, 'DefaultAxesFontSize', 14);
    set(0, 'DefaultAxesFontName', 'Arial');
    
end



function Drug_Intensity_percent_save = StatisticBoxPlot_v2(Drug_Intensity, Drug_PrecisionLow, Range_const1, Range_const2)

    if nargin <4
        Range_const2 = 2;
        if nargin < 3
            Range_const1 = 2;
            if nargin < 2
                Drug_PrecisionLow = 0;
            end
        end
    end
        
    Length_vec = length(Drug_Intensity);
    Length_real = min(sum(~isnan(Drug_Intensity)),sum(~isnan(Drug_Intensity)));
    
    % no idea why some cells has nan in one but not another channel.
    NotNaNIndex = find(~isnan(Drug_Intensity).*~isnan(Drug_Intensity));
    Drug_Intensity = Drug_Intensity(NotNaNIndex);
 
    % clean up NaN value when dealing with multiple kinds of objects
    
%% Statistics

%% Primary mask for rejecting extreme case
    Drug_raw  = Drug_Intensity;
    Drug_LowLim0 = -1;
    Drug_UpLim0 = 10;
    RatioMask_0 = (Drug_raw > Drug_LowLim0).*(Drug_raw < Drug_UpLim0);
    RatioMask_0 = RatioMask_0.*(Drug_Intensity > Drug_PrecisionLow);
     
    Drug_Intensity_0  = Drug_Intensity.*RatioMask_0;
    Drug_Intensity_0 = Drug_Intensity_0(Drug_Intensity_0~=0);
    
    Mean_Drug_Intensity_1 = mean(Drug_Intensity_0);
    Std_Drug_Intensity_1 = std(Drug_Intensity_0);

    
    %% Result masks 1 for abnormal results
%     Range_const1 = 2;
    
    % CD mask
    Drug_Intensity_LowLim1 = Mean_Drug_Intensity_1 - Range_const1*Std_Drug_Intensity_1;
    Drug_Intensity_UpLim1 = Mean_Drug_Intensity_1 + Range_const1*Std_Drug_Intensity_1;
    Drug_Intensity_mask1 = (Drug_Intensity > Drug_Intensity_LowLim1).*(Drug_Intensity < Drug_Intensity_UpLim1);
    
    %% Restatistics after picked out outliers
    Drug_Intensity_normal  = Drug_Intensity.*Drug_Intensity_mask1;
    Drug_Intensity_normal = Drug_Intensity_normal(Drug_Intensity_normal~=0);
   
    
    Mean_Drug_Intensity2 = mean(Drug_Intensity_normal);
    Std_Drug_Intensity2 = std(Drug_Intensity_normal);
    % Mean_off2 = mean(off_normal);
    % Std_off2 = std(off_normal);
 
    
    
    %% Result masks 2 for outliers (did not precisely measured due to the sample shape/movement)
%     Range_const2 = 2;   % pure culture
    
    % CD mask
    Drug_Intensity_LowLim2 = Mean_Drug_Intensity2 - Range_const2*Std_Drug_Intensity2;
    Drug_Intensity_UpLim2 = Mean_Drug_Intensity2 + Range_const2*Std_Drug_Intensity2;
    Drug_Intensity_mask2 = (Drug_Intensity > Drug_Intensity_LowLim2).*(Drug_Intensity < Drug_Intensity_UpLim2);
    
    % % off mask
    % off_LowLim2 = Mean_off2 - Range_const2*Std_off2;
    % off_UpLim2 = Mean_off2 + Range_const2*Std_off2;
    % off_mask2 = (off > off_LowLim2).*(off < off_UpLim2);
    
    
    %% Masked CD ratio result
%     CD_percent_raw  = (CD)./(CD+CH)*100;
    Drug_Intensity_percent_0 = Drug_raw.*Drug_Intensity_mask1;
    Drug_Intensity_percent_0 = Drug_Intensity_percent_0(Drug_Intensity_percent_0~=0);
    Drug_Intensity_percent_masked = Drug_raw.*Drug_Intensity_mask2.*(Drug_raw > Drug_LowLim0);
    Drug_Intensity_percent = Drug_Intensity_percent_masked(Drug_Intensity_percent_masked~=0);
    Drug_Intensity_percent_save = Drug_Intensity_percent_masked;
    Drug_Intensity_percent_save(Drug_Intensity_percent_save==0)= NaN;
    Drug_Intensity_percent_save = [Drug_Intensity_percent_save;NaN(Length_vec-Length_real,1)];
    
    %% Quick boxplot check
    figure;
    subplot(1,3,1);
    boxplot([Drug_raw]);
    hold on;
    scatter(ones(size(Drug_raw)).*(1+(rand(size(Drug_raw))-0.5)/10),Drug_raw,'r','filled')
    hold off;
%     ylim([-7 50]); 
    
    subplot(1,3,2);
    boxplot([Drug_Intensity_percent_0]);
    hold on;
    scatter(ones(size(Drug_Intensity_percent_0)).*(1+(rand(size(Drug_Intensity_percent_0))-0.5)/10),Drug_Intensity_percent_0,'r','filled')
    hold off;
    
    subplot(1,3,3);
    boxplot([Drug_Intensity_percent]);
    hold on;
    scatter(ones(size(Drug_Intensity_percent)).*(1+(rand(size(Drug_Intensity_percent))-0.5)/10),Drug_Intensity_percent,'r','filled')
    hold off;
    
    % print for quick quality check
    fprintf('\n');
    fprintf('CD percent std: %.3f \n',std(Drug_Intensity_percent));
    fprintf('CD percent mean: %.3f \n',mean(Drug_Intensity_percent));
    fprintf('Number of cell counted: %d \n', length(Drug_Intensity_percent));
    fprintf('Total cell number: %d \n\n', length(Drug_raw));
end