% Run in PositiveNegativeOverlays and SNPCalling folder
homedir = pwd;
SNPs_all = cell(0);
load('Coronavirus Targets.mat');
array_size = [sqrt(length(tar_names)), sqrt(length(tar_names))];

cd('../SampleData');
list = dir();
for n = 1:length(list)-2
    cd(list(n+2).name);

    comparison_mask = zeros(array_size);

    comparison_mask(3,2:5) = 1;
    comparison_mask(4,1:2) = 2;
    comparison_mask(4,3:5) = 3;
    comparison_mask(5,6) = 3;
    comparison_mask(5,1:3) = 4;
    comparison_mask(5,4:5) = 5;

    if ~exist('DATA.mat')
        cd(homedir);
        cd('../SampleData');
        continue;
    end
	load('DATA.mat');
    d = DATA{4};
    deriv_filter = zeros(size(d,1), size(d,2));

    for i = 3:size(deriv_filter,1)-2
        deriv_filter(i,:) = (1/12)*(-d(i+2,:)+8*d(i+1,:) - 8*d(i-1,:) + d(i-2,:));
    end

    deriv_max = max(deriv_filter(8:end,:));
    deriv_max = reshape(deriv_max, array_size).';

    SNPcalling = cell(0);
    for i = 1:max(max(comparison_mask))
        SNP_slopes = zeros(array_size);
        for j = 1:array_size(1)
            for k = 1:array_size(2)
                if comparison_mask(j,k) == i
                    SNP_slopes(j,k) = deriv_max(j,k);
                end
            end
        end

        [ind1,ind2] = find(SNP_slopes == max(max(SNP_slopes)));
        SNPcalling{end+1} = tar_names{array_size(2)*(ind1-1) + ind2};
    end
    fid = fopen('SNP Calls.txt', 'w');
	fprintf(fid, '%s\n', list(n+2).name);
    fprintf(fid, '%s\n', strjoin(SNPcalling,' '));
    SNPs_all{end+1} = SNPcalling;
    cd(homedir);
end

cd(homedir);
fid = fopen('SNP Calls_all.txt', 'w');
for n = 1:length(list)-2
    disp(list(n+2).name);
	fprintf(fid, '%s\n', list(n+2).name);
    if ~isempty(SNPs_all)
        disp(SNPs_all{n});
        fprintf(fid, '%s\n', strjoin(SNPs_all{n},' '));
    end
end
fclose(fid);