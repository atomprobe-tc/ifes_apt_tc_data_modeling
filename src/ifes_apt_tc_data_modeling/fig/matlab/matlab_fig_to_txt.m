% Markus Kühbach, 2023/04/28
% this converter is envisioned to be developed further to a generic converter
% which can transcode interactively made ranging definitions from Matlab
% figures to a text file to be parsed with e.g. the ifes_apt_tc_data_modeling/fig reader
% for ranging definitions in the research field of atom probe
clear;
format long;

fau_matlab_fig_range_file_names = { ...
    '02763-v01.RRNG.fig';
    'R04_11350-v01.rng.fig';
    'R04_15208-v04-clip-2_rng.fig';
    'R31_06413-v04-roi_rng.fig';
    'R56_01769.rng.fig';
    'R56_02476-v03_rrng.fig'};

for f = 1:length(fau_matlab_fig_range_file_names)
    disp(['Start converting ' fau_matlab_fig_range_file_names{f} '...']);

    % query the data from the figure
    tmp = openfig(fau_matlab_fig_range_file_names{f});
    fig = gcf;
    % assume that all relevant ranging definitions are represented as area
    % figure elements which have to be children of specific figure
    % components
    rng = findall(fig, 'Type', 'Area');
    iontype_id = 1;
    file_name = strcat(fau_matlab_fig_range_file_names{f}, '.txt');
    file_identifier = fopen(file_name, 'w');
    format_spec = '%s %d %d\n';
    for i = 1:1:length(rng)
        if rng(i).Type == 'area'
            ionname = rng(i).DisplayName;
            if ~strcmp(ionname, 'mass spectrum')
                disp(ionname)
                range = [rng(i).XData(1), rng(i).XData(length(rng(i).XData))];
              
                fprintf(file_identifier, format_spec, ionname, range(1), range(2));
                iontype_id = iontype_id + 1;
            end
        end
    end
    fclose(file_identifier);
    close(fig);
    disp(['Conversion ' fau_matlab_fig_range_file_names{f} ' success']);
end
disp('Conversion workflow completed');