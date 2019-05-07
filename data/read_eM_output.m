function [cMu, hL, nS, output_eM]=read_eM_output(filename)

%% Initialize variables.
delimiter = ':';
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(filename,'r');
output_text = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);

%% Create output variable
output_eM = table;
output_eM.Fields = cellstr(output_text{:, 1});
output_eM.Values = cellstr(output_text{:, 2});

%% Extract values of interest
cMu=str2num(output_eM.Values{match_str(output_eM.Fields,'Statistical Complexity')});
hL=str2num(output_eM.Values{match_str(output_eM.Fields,'History Length')});
nS=str2num(output_eM.Values{match_str(output_eM.Fields,'Number of Inferred States')});