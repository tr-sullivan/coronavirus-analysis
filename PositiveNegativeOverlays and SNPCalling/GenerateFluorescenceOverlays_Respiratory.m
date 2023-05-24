% Load in all info for experiments, targets, and panel layout
% Run in PositiveNegativeOverlays and SNPCalling folder
close all;
clear all;
homedir = pwd;
load('Respiratory Matrix.mat');
load('Coronavirus Targets.mat');

% Find all files to be analyzed
cd('../SampleData');
list = dir();

TIME = 0.5:0.5:30;
for i = 1:length(list)-2
    % Keep track of which dataset is being analyzed and move into its
    % directory
    name = list(i+2).name;
    name_split = split(name, '_');
    expt_comp = name_split{1};
    disp(['Experiment ', name]);
    
    % Account for multiple templates being present
    expected_probes = zeros(1,length(tar_names));
    target = split(expt_comp, ',');
    
    % Load in corresponding target matrix for each template
    for j = 1:length(target)
        id = find(strcmp(target{j}, expected(:,1)));
        for k = 1:length(expected_probes)
            new = expected{id,2}(k);
            old = expected_probes(k);
            if new == -1 || old == -1
                expected_probes(k) = -1;
            else
                expected_probes(k) = max(new, old);
            end
        end
    end
    
    % Robust vs. discriminating probes for coloring
    robustness = expected{1,2};
    
    % Load in analyzed data and calculate fluorescence derivative from
    % baseline-subtracted data
    cd(name);
    if ~exist('DATA.mat')
        cd(homedir);
        cd('../SampleData');
        continue;
    end
    load('DATA.mat');
    d_raw = DATA{1,2};
    d_norm = DATA{1,3};
    d_bsln = DATA{1,4};
    deriv_filter = zeros(size(d_bsln,1), size(d_bsln,2));
    for k = 3:size(deriv_filter,1)-2
        deriv_filter(k,:) = (1/12)*(-d_bsln(k+2,:)+8*d_bsln(k+1,:) - 8*d_bsln(k-1,:) + d_bsln(k-2,:));
    end
    for i = 3:size(deriv_filter,1)-2
        deriv_filter(i,:) = mean(deriv_filter((i-2):(i+2), :));
    end
    deriv_max = max(deriv_filter(8:end,:));

    % Filter expected amplifying probes based on choosing winners for each
    % SNP cluster
    amp_probes_id = find(expected_probes >= 1);
    amp_probes = tar_names(amp_probes_id);

    if any(strcmp(amp_probes, 'SARS-CoV-2'))
        sc2_snps = ["484";"417";"452";"501";"614"];
        for j = 1:length(sc2_snps)
            ids = find(contains(tar_names, sc2_snps{j}));
            for k = 1:length(ids)
                if ~ismember(ids(k), amp_probes_id)
                   expected_probes(ids(k)) = -2;
                end
            end
        end
    end

    % Initialize plotting
    handles = [];
    probes = cell(0);
    handles_ctrl = [];
    probes_ctrl = cell(0);
    handles_neg = [];
    probes_neg = cell(0);
    handles_negsnp = [];
    probes_negsnp = [];
    f = figure();
    for j = 1:length(tar_names)
        % If probe entry is -1, do not plot (i.e. for positive control spots)
        if expected_probes(j) == -1
            continue;
        % For probes expected to amplify, plot on the first set of axes and
        % color according to whether the probe of interest SNP-detecting
        % or not or if it is a control (extraction or amplification)
        elseif expected_probes(j) == 1
            subplot(1,2,1);
            if robustness(j) == 1
                handles = [handles, plot(TIME, d_bsln(:,j), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2)];
            elseif robustness(j) == 0
                handles = [handles, plot(TIME, d_bsln(:,j), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 2)];
            end
                probes{end+1} = tar_names{j};
            hold on;
        elseif expected_probes(j) == 2
            subplot(1,2,1);
            handles_ctrl = [handles_ctrl, plot(TIME, d_bsln(:,j), 'Color', 'g', 'LineStyle', '-', 'LineWidth', 2)];
            probes_ctrl{end+1} = tar_names{j};
            hold on;   
        % For probes not expected to amplify, plot on second set of axes.
        % Only color non-gray if the probe is required to be silent/weak as
        % part of the target set (e.g. for aac(6')-Ib)
        elseif expected_probes(j) == 0
            subplot(1,2,2);
            if robustness(j) == 0
                handles_neg = [handles_neg, plot(TIME, d_bsln(:,j), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 2)];
            elseif robustness(j) == 1
                handles_neg = [handles_neg, plot(TIME, d_bsln(:,j), 'Color', 'r', 'LineStyle', '-', 'LineWidth', 2)];
            elseif robustness(j) == -2
                handles_neg = [handles_neg, plot(TIME, d_bsln(:,j), 'Color', [0.7 0.7 0.7], 'LineStyle', '-', 'LineWidth', 2)];
            end
            probes_neg{end+1} = tar_names{j};
            hold on;
        elseif expected_probes(j) == -2
            subplot(1,2,2);
            handles_negsnp = [handles_negsnp, plot(TIME, d_bsln(:,j), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 2)];
            probes_negsnp{end+1} = tar_names{j};
            hold on;
        end
    end   
    for j = 1:length(probes_ctrl)
        probes{end+1} = probes_ctrl{1};
    end
    subplot(1,2,1);
    ys1 = ylim();
    xlabel('Time (min)');
    ylabel('Baseline-Subtracted Fluorescence Data');
    title('Expected Positive Probes');
    legend([handles, handles_ctrl], probes, 'Location', 'northwest');

    subplot(1,2,2);
    ys2 = ylim();
    if ys1(2) > ys2(2)
        ylim(ys1);
    else
        subplot(1,2,1);
        ylim(ys2);
        subplot(1,2,2);
    end
    xlabel('Time (min)');
    ylabel('Baseline-Subtracted Fluorescence Data');
    title('Expected Negative Probes');
    if ~isempty(probes_negsnp)
        legend(handles_negsnp, probes_negsnp, 'Location', 'northwest');
    end

    set(gcf,'PaperOrientation','landscape','WindowState','minimized');
    set(gcf, 'Position', get(0, 'Screensize'),'WindowState','minimized');
    print(['Exp ', name, '.pdf'], '-dpdf', '-fillpage')
    print(['Exp ', name, '.png'], '-dpng')
    close(f);

    cd(homedir);
    cd('../SampleData');
end

fclose('all');