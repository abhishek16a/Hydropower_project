function part2_save_single_figs_clear()
% Creates SINGLE, CLEAR, COMPARISON-WISE FIGURES for Part 2
% - For each input: panel L=Original MF, R=GA MF (no overlap)
% - For each output: panel L=Original MF, R=GA MF
% - For each control surface: panel L=Original, R=GA with matched color scale
%
% Files saved to ./part2_single_figs_clear

clc;
outdir = fullfile(pwd,'part2_single_figs_clear');
if ~exist(outdir,'dir'), mkdir(outdir); end

% --- Load baseline & optimised FIS
fis0  = readfis('AssistiveHomeFLC_Ext.fis');
fisGA = readfis('AssistiveHomeFLC_Ext_GA.fis');

% ---------------------- 1) INPUT MFs: side-by-side ----------------------
inNames = string({fis0.Inputs.Name});
for i = 1:numel(fis0.Inputs)
    fig = figure('Color','w','Position',[150 150 960 420]);
    t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile; hold on; grid on;
    plotmf(fis0,'input',i);
    title(sprintf('Input: %s (Original)', inNames(i)),'FontWeight','bold');
    xlabel('Universe'); ylabel('\mu'); legend('show','Location','best');

    nexttile; hold on; grid on;
    plotmf(fisGA,'input',i);
    title(sprintf('Input: %s (GA Optimised)', inNames(i)),'FontWeight','bold');
    xlabel('Universe'); ylabel('\mu'); legend('show','Location','best');

    title(t, sprintf('Input MF Comparison — %s', inNames(i)), 'FontWeight','bold');
    saveFig(fig, fullfile(outdir, sprintf('Input_%02d_%s_MF_compare.png', i, sanitize(inNames(i)))));
end

% ---------------------- 2) OUTPUT MFs: side-by-side ---------------------
outNames = string({fis0.Outputs.Name});
for i = 1:numel(fis0.Outputs)
    fig = figure('Color','w','Position',[150 150 960 420]);
    t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile; hold on; grid on;
    plotmf(fis0,'output',i);
    title(sprintf('Output: %s (Original)', outNames(i)),'FontWeight','bold');
    xlabel('Universe'); ylabel('\mu'); legend('show','Location','best');

    nexttile; hold on; grid on;
    plotmf(fisGA,'output',i);
    title(sprintf('Output: %s (GA Optimised)', outNames(i)),'FontWeight','bold');
    xlabel('Universe'); ylabel('\mu'); legend('show','Location','best');

    title(t, sprintf('Output MF Comparison — %s', outNames(i)), 'FontWeight','bold');
    saveFig(fig, fullfile(outdir, sprintf('Output_%02d_%s_MF_compare.png', i, sanitize(outNames(i)))));
end

% --------- 3) CONTROL SURFACES with matched color scale (zlim) ----------
% Define pairs: [inputA inputB -> outputIdx]
pairs = [1 3 1;    % Temp & Activity -> HeaterPower
         2 3 2;    % Light & Activity -> LightBrightness
         2 3 3;    % Light & Activity -> BlindsPosition
         5 4 4];   % CO2  & Humidity  -> VentilationFan
labels = ["Temp & Activity → HeaterPower", ...
          "Light & Activity → LightBrightness", ...
          "Light & Activity → BlindsPosition", ...
          "CO₂ & Humidity → VentilationFan"];

for k=1:size(pairs,1)
    inA = pairs(k,1); inB = pairs(k,2); outIdx = pairs(k,3);

    % Sample surfaces numerically to enforce identical clim/zlim
    [X0,Y0,Z0]   = local_gensurf_data(fis0,  inA, inB, outIdx);
    [X1,Y1,Z1]   = local_gensurf_data(fisGA, inA, inB, outIdx);
    zmin = min([Z0(:); Z1(:)]); zmax = max([Z0(:); Z1(:)]);

    fig = figure('Color','w','Position',[120 120 1100 460]);
    t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    surf(X0,Y0,Z0,'EdgeColor','none'); view(135,30); grid on; box on;
    title('Original'); xlabel(fis0.Inputs(inA).Name); ylabel(fis0.Inputs(inB).Name);
    zlabel(fis0.Outputs(outIdx).Name);
    colormap(gca,'parula'); caxis([zmin zmax]); colorbar;

    nexttile;
    surf(X1,Y1,Z1,'EdgeColor','none'); view(135,30); grid on; box on;
    title('GA Optimised'); xlabel(fisGA.Inputs(inA).Name); ylabel(fisGA.Inputs(inB).Name);
    zlabel(fisGA.Outputs(outIdx).Name);
    colormap(gca,'parula'); caxis([zmin zmax]); colorbar;

    title(t, labels(k), 'FontWeight','bold');
    saveFig(fig, fullfile(outdir, sprintf('Surface_%02d_compare.png', k)));
end

fprintf('Saved clear comparison figures in: %s\n', outdir);

% ------------------------- helpers -------------------------
function s = sanitize(str)
    s = regexprep(char(str),'[^A-Za-z0-9_-]','_');
end

function saveFig(fig, filename)
    try
        exportgraphics(fig, filename, 'Resolution', 300);
    catch
        saveas(fig, filename);
    end
    close(fig);
end

function [X,Y,Z] = local_gensurf_data(fis, inA, inB, outIdx)
    % Build grids from actual variable ranges
    rA = fis.Inputs(inA).Range; rB = fis.Inputs(inB).Range;
    X  = linspace(rA(1), rA(2), 41);
    Y  = linspace(rB(1), rB(2), 41);
    [XX,YY] = meshgrid(X,Y);
    % build samples, evaluate only once, then reshape
    n = numel(XX);
    % make a default input row per sample (fill non-used inputs with midpoints)
    m = numel(fis.Inputs);
    U = zeros(n,m);
    for j=1:m
        rr = fis.Inputs(j).Range;
        U(:,j) = (rr(1)+rr(2))/2; % midpoint
    end
    U(:,inA) = XX(:);
    U(:,inB) = YY(:);
    YYhat = evalfis(fis,U);
    Z = reshape(YYhat(:,outIdx), size(XX));
end
end
