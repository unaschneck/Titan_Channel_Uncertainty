clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE SINUOSITY OF TITAN RIVERS (in chunks) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath(fullfile('functions')))                                    % all functions for extract and analysis channel data
addpath(genpath('data'))                                                   % raw flyby datas and marked feature data

% where is the data with the shape files
Titan_rivers_dir = fullfile('data','TITAN_MAPPED_CHANNELS');                          % location of directory with all Titan mapped river shape files
Titan_river = dir(Titan_rivers_dir);                                                  % all sub folders containing specific river data

addpath(genpath(Titan_rivers_dir));

fprintf('There are %i river features in %s\n',numel(Titan_river)-2,Titan_rivers_dir); % minus 2 to ignore folders '.' and '..'

results_dir = 'results';
if ~isfolder(results_dir)
    mkdir(results_dir);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('==================================================================')
disp('CALCULATING SINUOISTY')
disp('==================================================================')

max_subdivisions = 20;                                                     % maximum number of subdivision within the channel to consider
k_sin = 1;                                                                 % counter for sinuosity calculation in inner loop

figure('Name','Channel sinuosity');
hold on;
grid on
xlabel('decreasing number of subdivisions')
ylabel('sinusoity')


line_color = distinguishable_colors(numel(Titan_river));
most_sin = 1;
windy_river = '';

for r = 1:numel(Titan_river)   
    if Titan_river(r).isdir && ~ismember(Titan_river(r).name , {'.','..'}) % ignore folders '.' and '..'

        path2data = fullfile(Titan_rivers_dir, Titan_river(r).name);
        all_mat = dir(fullfile(path2data, '*.mat'));

        fprintf('Feature: %s\n',Titan_river(r).name)

        if ~isempty(all_mat)

            for j = 1:length(all_mat)

                fname = all_mat(j);
                fprintf('Calculating sinuosity for: %s\n', fname.name);

                avg_sin(k_sin,:) = calc_average_sinuosity(fname.name,max_subdivisions);

                plot(avg_sin(k_sin,:),'LineWidth',2,'Color',line_color(r,:)) % the same river but from different flybys will be the same color (this always returns the same sinuosity, is that reasonable? -- geomtry structure containing the centerline is the same for same feature across different flybys, does that make sense or should it be slightly different?)

                if max(avg_sin(k_sin,:)) > most_sin
                    most_sin = max(avg_sin(k_sin,:));
                    windy_river = fname.name;
                end
                k_sin = k_sin + 1;

            end
        else
            fprintf('no .mat file available for %s\n',Titan_river(r).name)
        end
    end
end
hold off;

fprintf('The most sin river is %s with a maximum sin of %0.2f\n',windy_river,most_sin)

avg_sin_global = avg_sin(:,end); % using the entire length for straight line distance

% want the bin edges to be the same for the different groups in the histogram
bin_edges = linspace(min([avg_sin(:); avg_sin_global]), max([avg_sin(:); avg_sin_global]), 30);

figure('Name','Channel sinuosity distribution');
histogram(avg_sin,bin_edges,'DisplayName','All segments')
hold on
histogram(avg_sin_global,bin_edges,'DisplayName','One segment')
xlabel('sinuosity')
ylabel('count')
legend('show')
grid on;
title(sprintf('Total number of %i maps for %i unique features',k_sin,numel(Titan_river)-2))



