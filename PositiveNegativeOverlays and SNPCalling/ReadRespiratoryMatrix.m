% Reads an excel file, respiratory matrix.xlsx, that specifies what probes are
% expected (1), not expected (0), or not to be plotted (-1). Curves that
% are not plotted include positive controls or irrelevant discriminatory
% probes (as these could be amplifying or not but the output is not
% important for differentiation of this particular strain)

fname = 'multiple pathogen matrix.xlsx';
sheets = sheetnames(fname);

expected = cell(length(sheets), 2);
for i = 1:length(sheets)
    expected{i,1} = char(sheets(i));
    matrix = xlsread(fname, sheets(i));
    matrix = reshape(matrix', 1, 36);
    expected{i,2} = matrix;
end

save('Respiratory Matrix.mat', 'expected');