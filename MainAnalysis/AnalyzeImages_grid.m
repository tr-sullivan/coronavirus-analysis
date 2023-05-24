close all
clear all
%% Set parameters
radius_range = [10 20]; %range of acceptable radii for spot finding
rnd = -1; %used for rounding
int_corr = 0.2; %used for re-scaling image for masking
border = 20; %number of pixels from each edge for zoomed-in crop of array
threshold = 15; %15 percent of the range used to detect T_t
threshold_deriv = 50; %threshold for detecting amplification (in terms of fluorescence derivative)
threshold_delta = 10; %threshold for detecting amplification (in terms of total change in fluorescence)
threshold_outlier = 4; %number of consecutive time points of the derivative at least 2 standard deviations from the mean

%load spotfile with probe IDs  
load('../../MainAnalysis/Coronavirus Targets.mat');
array_size = [sqrt(length(tar_names)), sqrt(length(tar_names))]; %[number of rows, number of columns]
shifts = zeros(sqrt(length(tar_names))-1,1);

%set mask for control and non-control spots
target_mask = ones(1, length(tar_names));
target_mask_NC = zeros(1, length(tar_names));
target_mask_PC = zeros(1, length(tar_names));
for i = 1:length(tar_names)
    if strcmp(tar_names{i}, 'NC')
        target_mask(i) = 0;
        target_mask_NC(i) = 1;
    end
    if strcmp(tar_names{i}, 'PC')
        target_mask(i) = 0;
        target_mask_PC(i) = 1;
    end
end
%% Perform file extraction and deletion if necessary
% get all raw files in current working directory
list=dir('#*.raw');

num_img = floor(length(list));
TIME = [0.5:0.5:num_img/2];

IMAGES = cell(num_img,1);
for i = 1:num_img
    % read in data
    fid = fopen(list(i).name, 'r');
    if fid == -1
        error('Cannot open file: %s', FileName);
    end
    data = fread(fid, [780,780], 'uint16');
    fclose(fid);
    data = reshape(data, [780, 780]);

    % convert to intensity
    data(data < 0) = 0;
    I = mat2gray(data);

    % write adjusted images to PNG
    Im_out = imadjust(I, [0 0.5], [0 1]);
    fname = ['S', sprintf('%02d', str2num(extractBetween(string(list(i).name), '#', '_'))),'_', char(extractBetween(list(i).name,'_', '.')),'.jpg'];
    imwrite(Im_out, fname);

    id = str2num(extractBetween(string(list(i).name), '#', '_'));
    IMAGES{id,1}=data;
end
%write subtracted image
sub_2 = imread('S02_780x780.jpg');
sub_end = imread(['S', num2str(length(TIME)), '_780x780.jpg']);
sub = sub_end - sub_2;
sub_out = imadjust(sub, [0 0.2], [0 1]);
colormap gray;
imagesc(sub_out);
print('Sfinal_S2', '-dpdf', '-bestfit');
save('Images.mat','IMAGES');

%% Masking 
%make array's masks
if ~exist('Masks.mat', 'file')
    new_masks = 1;
else
    new_masks = menu('overwrite previous masks?', 'yes','no');
end
if new_masks == 1
    rgb1 = imadjust(mat2gray(IMAGES{end}), [0 0.5], [0 1]); 
    [MASKS(:,1), MASKS_BG(:,1)] = find_masks_grid(rgb1, int_corr, radius_range, rnd, array_size, shifts);
    [target_number,~] = size(MASKS(:,1));
    save('Masks.mat','MASKS', 'MASKS_BG');
    save('target_number.mat','target_number');
    close all
else 
    load('target_number.mat');
    load('Masks.mat');
end

%crop image to +/- border pixels from the edge of the array and save
mask = zeros(size(MASKS{1}));
for j = 1:length(MASKS)
    mask = mask + MASKS{j};
end
[r,c] = find(mask);
r1 = min(r);
r2 = max(r);
c1 = min(c);
c2 = max(c);
sub = 5*(IMAGES{end} - IMAGES{2});
sub_roi = sub(max(r1-border, 1):min(r2+border, size(IMAGES{1},1)), max(c1-border,1):min(c2+border, size(IMAGES{1},2)));
image_path = [pwd,'\ArrayCrop','.png'];
imwrite(imadjust(uint16(sub_roi), [0 0.3], [0 1]), image_path);

fprintf('wait...\n');

%calculate each spot's raw intensity
DATA = cell(1);
DATABG = cell(1);
  
%target intetsity
for i = 1:num_img
    disp(i);
    for ii = 1:target_number
        IMAGE = IMAGES{i};
        MASK = MASKS{ii};
        MASKBG = MASKS_BG{ii};
        DATA{1,1}(i,ii) = sum(sum(double(IMAGE) .* double(MASK))) / sum(double(MASK(:)));
        DATABG{1,1}(i,ii) = sum(sum(double(IMAGE) .* double(MASKBG))) / sum(double(MASKBG(:)));
        DATA{1,2}(i,ii) = DATA{1,1}(i,ii) - DATABG{1,1}(i,ii);
    end
end

s = size(DATA{1,1});

save ('raw_DATA.mat','DATA','DATABG');
save ('TIME.mat','TIME');

%% Normalization
load ('raw_DATA.mat');
load ('target_number.mat');

[a, b] = size (DATA{1,1}); 

%find maximum of all non-control data points
max_val = max(max((DATA{1,2}(2:end,:) - mean(DATA{1,2}(3:7,:))).* target_mask));

% find maximum and minimum of all control data points during baseline
ctrl_PC = DATA{1,2} .* target_mask_PC;
ctrl_NC = DATA{1,2} .* target_mask_NC;
ctrl_PC_baseline = ctrl_PC(3:7,:);
ctrl_NC_baseline = ctrl_NC(3:7,:);
max_val_ctrl = median(max(ctrl_PC_baseline(ctrl_PC_baseline ~= 0)));
min_val_ctrl = median(min(ctrl_NC_baseline(ctrl_NC_baseline ~= 0)));

if isnan(max_val_ctrl)
    max_val_ctrl = median(max(ctrl_NC_baseline(ctrl_NC_baseline ~= 0)));
end
if isnan(min_val_ctrl)
    min_val_ctrl = median(min(ctrl_PC_baseline(ctrl_PC_baseline ~= 0)));
end

%normalization to 0-1
for k = 1:target_number
    if strcmp(tar_names{k,1}, 'PC') || strcmp(tar_names{k,1}, 'NC')
        DATA{1,3}(:,k) = normalization_current(DATA{1,2}(:,k), tar_names{k}, [min_val_ctrl, max_val_ctrl]);
    else
        DATA{1,3}(:,k) = normalization_current(DATA{1,2}(:,k), tar_names{k}, max_val);
    end   
end
DATA{1,4} = DATA{1,2} - mean(DATA{1,2}(3:7,:));

save ('DATA.mat','DATA');
save ('TIME.mat','TIME');

fprintf('Done.\n');

%% Amplification Detection
d = DATA{4};
deriv_filter = zeros(size(d,1), size(d,2));
for i = 3:size(deriv_filter,1)-2
    deriv_filter(i,:) = (1/12)*(-d(i+2,:)+8*d(i+1,:) - 8*d(i-1,:) + d(i-2,:));
end
deriv_baseline = mean(deriv_filter(5:10,:));
deriv_sd = std(deriv_filter(5:10,:));

longestchain = zeros(1,target_number);
for i = 1:target_number
    bl = deriv_baseline(i);
    sd = deriv_sd(i);
    dat = deriv_filter(:,i);
    outliers = dat > bl + 1.5*sd;
    count = 0;
    countmax = 0;
    for j = 1:length(outliers)
        if outliers(j)
            count = count + 1;
        else
            countmax = max([countmax, count]);
            count = 0;
        end
    end
    countmax = max([countmax, count]);
    longestchain(i) = countmax;
end
deriv_max = max(deriv_filter(8:end,:));
deriv_max = reshape(deriv_max, array_size).';
delta_RFU = d(end,:) - mean(d(3:7,:));
delta_RFU = reshape(delta_RFU, array_size).';
% amp = (deriv_max > threshold_deriv) & (delta_RFU > threshold_delta);
amp = reshape(longestchain >= threshold_outlier, [array_size])' & delta_RFU > threshold_delta;

plotform = amp';
plotform = plotform(:);

%% Plot Raw Data
% subplot the data
nCol = array_size(2);
nRow = array_size(1);
Rect = [0.07, 0.07, 0.90, 0.88];
space_x = 0.0001;
space_y = 0.001;
AxisPos = moPlotPos(nCol, nRow, Rect, space_x, space_y);
set(gcf, 'Position', [100, 100, 600, 600*3/3.5],'WindowState','minimized');

plot_pos = [100,500;...
            500,500;...
            100,50;...
            500,50];

max_raw = max(max(DATA{1,1}));
min_raw = min(min(DATA{1,1}));
figure();
set(gcf, 'Position', [plot_pos(1,:), 600, 600*3/3.5],'WindowState','minimized');

for i = 1:target_number
    Data_ct = DATA{1,2}(:,i);
    axes('Position', AxisPos(i, :));
    if plotform(i) == 0
        color = 'c';
    elseif plotform(i) == 1
        color = 'r';
    else
        color = '-.r';
    end
    plot(TIME, Data_ct, color, 'linewidth', 2); 
    hold on
    set(gca, 'Color', 'None');
    y_max = max(Data_ct);
    ylim([0, max(max(DATA{1,2}))]);
    %ylim([0, 5000]);

    ys = ylim();
    text (3,ys(1) + (ys(2) - ys(1))*0.8,tar_names{i,1}, 'FontSize', 8);

    set(gca, 'xlim',([TIME(1)-0.5 TIME(end)]));
    set(gca, 'XTick',[0:10:TIME(end)])

    % plot x values only on specific subplots 
    x_array = [b-nCol+1:b];
    if ismember(i,x_array)
        set(gca,'XTickLabel',[0:10:TIME(end)],'FontSize',8);
        xlabel('Time (min)','FontSize',8)
    else
        set(gca,'XTickLabel',[]);
    end

    set(gca,'xcolor','k');
    set(gca,'ycolor','k');
end

set(gcf,'PaperOrientation','landscape');
set(gcf, 'Position', get(0, 'Screensize'));
print('Results_raw_Amp', '-dpng')


%% Plot Baseline-Subtracted Data
% subplot the data
nCol = array_size(2);
nRow = array_size(1);
Rect = [0.07, 0.07, 0.90, 0.88];
space_x = 0.0001;
space_y = 0.001;
AxisPos = moPlotPos(nCol, nRow, Rect, space_x, space_y);
% set(gcf, 'Position', [100, 100, 600, 600*3/3.5],'WindowState','minimized');
set(gcf, 'Position', [500, 500, 400, 400*3/3.5],'WindowState','minimized');

plot_pos = [100,500;...
            500,500;...
            100,50;..., 
            500,50];

max_raw = max(max(DATA{1,1}));
min_raw = min(min(DATA{1,1}));
figure();
set(gcf, 'Position', [plot_pos(1,:), 600, 600*3/3.5],'WindowState','minimized')

for i = 1:target_number
    Data_ct = DATA{1,4}(:,i);
    axes('Position', AxisPos(i, :));
    if plotform(i) == 0
        color = 'c';
    elseif plotform(i) == 1
        color = 'r';
    else
        color = '-.r';
    end
    plot(TIME, Data_ct, color, 'linewidth', 2); 
    hold on
    set(gca, 'Color', 'None');
    y_max = max(Data_ct);
    ylim([min(min(DATA{1,4})), max(max(DATA{1,4}))]);

    ys = ylim();
    text (3,ys(1) + (ys(2) - ys(1))*0.8,tar_names{i,1}, 'FontSize', 8);

    set(gca, 'xlim',([TIME(1)-0.5 TIME(end)]));
    set(gca, 'XTick',[0:10:TIME(end)])

    % plot x values only on specific subplots 
    x_array = [b-nCol+1:b];
    if ismember(i,x_array)
        set(gca,'XTickLabel',[0:10:TIME(end)],'FontSize',8);
        xlabel('Time (min)','FontSize',8)
    else
        set(gca,'XTickLabel',[]);
    end

    set(gca,'xcolor','k');
    set(gca,'ycolor','k');
end

set(gcf,'PaperOrientation','landscape');
set(gcf, 'Position', get(0, 'Screensize'));
print('Results_baseline_subtract_Amp', '-dpng')

%% Plot Normalized Data
figure();
set(gcf, 'Position', [plot_pos(1,:), 600, 600*3/3.5],'WindowState','minimized');

for i = 1:target_number
    Data_ct = DATA{1,3}(:,i);
    axes('Position', AxisPos(i, :));
    if plotform(i) == 0
        color = 'c';
    elseif plotform(i) == 1
        color = 'r';
    else
        color = '-.r';
    end
    plot(TIME, Data_ct, color, 'linewidth', 2); 
    hold on
    set(gca, 'Color', 'None');
    y_max = max(Data_ct);
    ylim([min(0, min(Data_ct)), max(1, max(Data_ct)) + 0.1]);

    [Tt(i), threshold_val(i)] = find_Tt_current( Data_ct, TIME*60, y_max, threshold, tar_names{i});
    threshold_times(i,1) = Tt(i);

    [Te(i), threshold_upper(i)] = find_Te_current(Data_ct, TIME*60, y_max, threshold, tar_names{i});
    threshold_times(i,2) = Te(i);

    ys = ylim();
    if ~isnan(Tt(i)) && plotform(i)>0
        line ([1 TIME(end)], [threshold_val(i) threshold_val(i)],'Color',[0.2 0.2 0.2]);
        text (3,ys(1) + (ys(2) - ys(1))*0.2,['Tt = ' num2str(Tt(i))], 'FontSize', 8);
    end
    text (3,ys(1) + (ys(2) - ys(1))*0.8,tar_names{i,1}, 'FontSize', 8);
 
    set(gca, 'xlim',([TIME(1)-0.5 TIME(end)]));
    set(gca, 'XTick',[0:10:TIME(end)]);

    % plot x values only on specific subplots 
    x_array = [b-nCol+1:b];
    if ismember(i,x_array)
        set(gca,'XTickLabel',[0:10:TIME(end)],'FontSize',8);
        xlabel('Time (min)','FontSize',8);
    else
        set(gca,'XTickLabel',[]);
    end

    y_array = [1:nCol:b];
    if ismember(i,y_array)
        set(gca,'YTickLabel',[0 0.5 1],'FontSize',8);
    else
        set(gca,'YTickLabel',[]);
    end

    set(gca,'xcolor','k');
    set(gca,'ycolor','k');
end

chip_name = ['chip-1 qPCR Data'];
saveas(gcf, chip_name, 'pdf')

set(gcf,'PaperOrientation','landscape','WindowState','minimized');
set(gcf, 'Position', get(0, 'Screensize'),'WindowState','minimized');
print('Results_multi_bg_Amp', '-dpdf', '-fillpage')
print('Results_norm_Amp', '-dpng')

xlswrite('Data.xls', DATA{1,1}, 1);
xlswrite('Data_norm.xls', DATA{1,3}, 1);
Tt = Tt';
writetable(table(tar_names, Tt), 'Thresholds.xlsx');
writetable(table(tar_names, longestchain'), 'DetectedProbes.xlsx');