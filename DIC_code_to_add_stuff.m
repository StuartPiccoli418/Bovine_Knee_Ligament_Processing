clear;
clc;
close all;

%% Set Input/Output Folders
%Make sure each file has 32 columns or it will be skipped
%No: Du/dt, Du/^2/d^2t, Filename,
%Make sure analog data is included with the data (0, Load, 2 etc)
%Format: x,y,z,u,v,w,d,x,y,z,u,v,w,d,count,time_0,Dev(1-7)/ai(1-7),0,load,2,3,4,5,6,7
%Note XX = X' YY = Y' ZZ = Z', Make sure the file is in .csv
%Load Constant = 1865.5

inputFolder = 'C:\Users\stuar\Downloads\DIC_Input'; %change to input root file path of choice
outputFolder = 'C:\Users\stuar\Downloads\DIC_output'; %Change to output root file path of choice

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Process All Files in Folder
fileList = dir(fullfile(inputFolder, '*.csv'));

summaryIndex = 1;
for k = 1:length(fileList)

    %% Read file
    inputFile = fullfile(inputFolder, fileList(k).name);
    data = readtable(inputFile, 'PreserveVariableNames', true, 'ReadVariableNames', true);
    disp(['📄 File: ', fileList(k).name]);

    %% Check column count
    if width(data) ~= 32
        warning('⚠️ File %s has %d columns, expected 32. Skipping.', fileList(k).name, width(data));
        continue;
    end

    %% Rename columns
    data.Properties.VariableNames = {
        'X_mm_1', 'Y_mm_1', 'Z_mm_1', 'U_mm_1', 'V_mm_1', 'W_mm_1', 'D_mm_1', ...
        'X_mm_2', 'Y_mm_2', 'Z_mm_2', 'U_mm_2', 'V_mm_2', 'W_mm_2', 'D_mm_2', ...
        'Count', 'Time_0', 'Dev1_ai0', 'Dev1_ai1', 'Dev1_ai2', 'Dev1_ai3', ...
        'Dev1_ai4', 'Dev1_ai5', 'Dev1_ai6', 'Dev1_ai7', 'Zero', ...
        'Load', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven'
        };

    %% Create Displacement Columns
    data.XX_mm = data.X_mm_1 + data.U_mm_1;
    data.YY_mm = data.Y_mm_1 + data.V_mm_1;
    data.ZZ_mm = data.Z_mm_1 + data.W_mm_1;

    data.XX2_mm = data.X_mm_2 + data.U_mm_2;
    data.YY2_mm = data.Y_mm_2 + data.V_mm_2;
    data.ZZ2_mm = data.Z_mm_2 + data.W_mm_2;

    %% Compute Magnitude and Strain
    data.Magnitude = sqrt((data.XX_mm - data.XX2_mm).^2 + ...
        (data.YY_mm - data.YY2_mm).^2 + ...
        (data.ZZ_mm - data.ZZ2_mm).^2);

    Magnitude_first = data.Magnitude(2);
    data.Strain = (data.Magnitude - Magnitude_first) / Magnitude_first;

    %% Warning for high strain
    if max(data.Strain) > 1
        fprintf('⚠️  High strain in %s: %.2f%%\n', fileList(k).name, max(data.Strain) * 100);
    else
        fprintf('✓  Strain OK in %s: %.2f%%\n', fileList(k).name, max(data.Strain) * 100);
    end

    %% Create Load_N
    data.Load_N = data.Load * 1865.5;

    %% Print Max Load_N
    fprintf('Max Load_N is %.2f\n', max(data.Load_N))

    %% Input Area 
    Area = input('Enter the correct area: ');
    data.Stress = data.Load_N ./ Area;

    %% Rearranging columns — simple and reliable way
    % Define the desired order (use names)
    columnOrder = {
        'X_mm_1', 'Y_mm_1', 'Z_mm_1', 'U_mm_1', 'V_mm_1', 'W_mm_1', ...
        'XX_mm', 'YY_mm', 'ZZ_mm', ...
        'X_mm_2', 'Y_mm_2', 'Z_mm_2', 'U_mm_2', 'V_mm_2', 'W_mm_2', ...
        'XX2_mm', 'YY2_mm', 'ZZ2_mm', ...
        'D_mm_1', 'D_mm_2', 'Count', 'Time_0', ...
        'Dev1_ai0', 'Dev1_ai1', 'Dev1_ai2', 'Dev1_ai3', ...
        'Dev1_ai4', 'Dev1_ai5', 'Dev1_ai6', 'Dev1_ai7', ...
        'Load', 'Load_N', 'Magnitude', 'Stress', 'Strain', ...
        'One', 'Two', 'Three', 'Four', 'Five', 'Six', 'Seven'
        };

    % Reorder if all names exist
    commonCols = intersect(columnOrder, data.Properties.VariableNames, 'stable');
    data = data(:, commonCols);  % Keeps columns that exist, in desired order

    %% Loop through files

    % Settings
    threshold = 1;
    limit = 0.05;
    outputFolder = 'C:\Users\stuar\Downloads\DIC_output';

    % Ensure output folder exists
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Load strain data (already loaded in your loop)
    y = data.Strain;
    x = 0.2 * (0:length(y)-1);
    z = data.Stress;

    % Force row vectors
    x = x(:)';
    y = y(:)';
    z = z(:)';

    % Skip if all values are below limit
    if all(y < limit)
        disp(['Skipping plot for "' fileList(k).name '" — yline is above all data']);
        continue;
    end
    %% Pick Two Points:
    % First point is around 0.05 Strain,
    % Second point is around peak of Strain
    % Graph will auto skip if all points are below 0.05

    % Load (N) vs Time (0.2 Seconds)
    fig1 = figure;
    plot(x, data.Load_N, '-.k.')
    xlabel('Time (0.2 Seconds)');
    ylabel('Force (N)');
    title(fileList(k).name, 'Load Force', 'Interpreter', 'none');
    grid on;
    hold on;
    xlims = xlim; ylims = ylim;
    x_text = xlims(1) + 0.04 * (xlims(2) - xlims(1));
    y_text = ylims(2) - 0.20 * (ylims(2) - ylims(1));
    text(x_text, y_text, sprintf('Max Load = %.2f N', max(data.Load_N)), ...
    'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold');
    close(fig1)

    % Stress v Strain Graph
    fig = figure;
    plot(y, z, '-.b.');
    xlabel('Strain');
    ylabel('Stress');
    title(fileList(k).name, 'Interpreter', 'none');
    grid on;
    hold on;
    xline(limit, 'k-', {'Toe Region', 'y=0.05'}, 'LineWidth', 1);

    % Let user manually select zoom region
    [xp, yp] = ginput(2);  % User selects zoom region

    % Only define ranges after valid input
    if numel(xp) ~= 2 || numel(yp) ~= 2 || ...
            any(isnan(xp)) || any(isnan(yp)) || ...
            (diff(yp) == 0) || (diff(xp) == 0)
        disp('⚠️ Invalid or skipped zoom selection — skipping this plot.');
        close(fig);
        continue;
    end

    x_range = sort(xp);
    y_range = sort(yp);

    % Now that ranges are defined, apply filters
    in_range = (x >= x_range(1) & x <= x_range(2)) & ...
        (y >= y_range(1) & y <= y_range(2));
    valid_points = in_range & y <= thre13shold;

    x_fit = x(valid_points);
    y_fit = y(valid_points);

    % Validate click
    if numel(yp) == 2 && all(isfinite(yp))
        x_range = sort(xp);
        y_range = sort(yp);  % Ensures increasing order
        if numel(y_range) == 2 && all(isfinite(y_range)) && diff(y_range) > 0
            ylim(y_range);
        else
            warning('⚠️ Invalid y_range: skipping ylim setting.');
            close(fig)
        end
        xlim(x_range);
    else
        warning('Invalid zoom range selected. Skipping.');
        close(fig)
        continue;
    end

    % Filter data points in selected region
    in_range = (x >= x_range(1) & x <= x_range(2)) & ...
        (y >= y_range(1) & y <= y_range(2));

    x_fit = x(in_range);
    y_fit = y(in_range);

    fprintf('Total points in selected region: %d\n', numel(x_fit))
    fprintf('Valid x_fit: %d\n', sum(~isnan(x_fit)))
    fprintf('Valid y_fit: %d\n', sum(~isnan(y_fit)))
     if numel(x_fit) < 2 || any(isnan(x_fit)) || any(isnan(y_fit))
        disp('⚠️ Not enough valid points for best-fit line. Skipping.');
        summary(summaryIndex).Filename = fileList(k).name;
        summary(summaryIndex).Status = 'Failed';
        summaryIndex = summaryIndex + 1;
        close(fig);
        continue;
     end


    % Fit linear regression
    coeffs = polyfit(x_fit, y_fit, 1);
    slope = coeffs(1);
    intercept = coeffs(2);
    disp(['Best-fit slope: ' num2str(coeffs(1))]);

    x_line = linspace(min(x_fit), max(x_fit), 100);
    y_line = polyval(coeffs, x_line);
    plot(x_line, y_line, 'r-', 'LineWidth', 2);

    % Update limits AFTER plotting everything
    xlims = xlim;
    ylims = ylim;

    % Slope equation
    range_x = xlims(2) - xlims(1);
    range_y = ylims(2) - ylims(1);
    x_text_eqn = xlims(1) + 0.04 * range_x;
    y_text_eqn = ylims(2) - 0.20 * range_y;
    eqnStr = sprintf('y = %.6fx + %.6f', slope, intercept);
    text(x_text_eqn, y_text_eqn, eqnStr, 'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold');

    % Max strain text
    maxStrain = max(data.Strain);
    x_text_max = xlims(1) + 0.04 * range_x;
    y_text_max = ylims(2) - 0.27 * range_y;
    text(x_text_max, y_text_max, sprintf('Max Strain = %.2f%%', maxStrain * 100), ...
        'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold');

    % === Set legend LAST ===
    legendEntries = {'Data'};
    if exist('x_line', 'var')
        legendEntries{end+1} = 'Best-Fit Line';
    end
    legend(legendEntries, 'Location', 'north');


    % Step 1: Fit the linear regression
    y_fit_pred = polyval(coeffs, x_fit);  % Predicted y values

    % Step 2: Calculate R² (coefficient of determination)
    SS_res = sum((y_fit - y_fit_pred).^2);                 % Residual sum of squares
    SS_tot = sum((y_fit - mean(y_fit)).^2);                % Total sum of squares
    R_squared = 1 - (SS_res / SS_tot);                     % R²

    % R² annotation
    x_text_r2 = xlims(1) + 0.05 * range_x;
    y_text_r2 = ylims(2) - 0.10 * range_y;
    infoStr = sprintf('Slope = %.5f\nR^2 = %.4f', coeffs(1), R_squared);
    text(x_text_r2, y_text_r2, infoStr, ...
        'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold', 'BackgroundColor', 'w');

    % Step 3: Display on the plot
    xlims = xlim;
    ylims = ylim;

    % Position text near top-left corner
    x_text = xlims(1) + 0.05 * (xlims(2) - xlims(1));
    y_text = ylims(2) - 0.10 * (ylims(2) - ylims(1));

    % Slope and R² text
    infoStr = sprintf('Slope = %.5f\nR^2 = %.4f', coeffs(1), R_squared);

    text(x_text, y_text, infoStr, ...
        'FontSize', 12, 'Color', 'k', 'FontWeight', 'bold', 'BackgroundColor', 'w');

    % Save Stress v Strain to output folder
    [~, nameOnly, ~] = fileparts(fileList(k).name);
    filename = fullfile(outputFolder, [nameOnly '_StressStrain.png']);
    exportgraphics(fig, filename, 'Resolution', 300);

    % Save Load(N) v Time to output folder
    filename1 = fullfile(outputFolder, [nameOnly '_LoadTime.png']);
    exportgraphics(fig1, filename1, 'Resolution', 300);
    close(fig1)
    
    % Save processed data
    newFileName = [nameOnly, '_Processed.csv'];
    outputFile = fullfile(outputFolder, newFileName);
    writetable(data, outputFile);

    % Write to Log File
    summary(summaryIndex).Filename = fileList(k).name;
    summary(summaryIndex).MaxStress = max(data.Stress);
    summary(summaryIndex).MaxStrain = max(data.Strain);
    summary(summaryIndex).Slope = coeffs(1);
    summary(summaryIndex).R_squared = R_squared;
    summary(summaryIndex).eqnStr = eqnStr;
    summary(summaryIndex).Load_N = max(data.Load_N);
    summary(summaryIndex).Area = Area;
    summary(summaryIndex).Status = 'Success';
    summaryIndex = summaryIndex + 1;
end

disp('✅ All files processed and saved.');
summaryTable = struct2table(summary);
writetable(summaryTable, fullfile(outputFolder, 'SummaryLog.csv'));


