function fis = part1_build_and_show_fis(saveAs, savePngs)
% PART 1 — Build and visualize an assistive home Mamdani FIS (no "Any" MF)
% Usage:
%   fis = part1_build_and_show_fis();                 % default save name, no PNGs
%   fis = part1_build_and_show_fis('myFIS.fis',true); % save figures as PNGs

if nargin<1 || isempty(saveAs),  saveAs = 'AssistiveHomeFLC_Ext.fis'; end
if nargin<2 || isempty(savePngs), savePngs = false; end

%% ---------- Build Mamdani FIS ----------
fis = mamfis('Name','AssistiveHomeFLC_Ext', ...
    'AndMethod','min','OrMethod','max', ...
    'ImplicationMethod','min','AggregationMethod','max', ...
    'DefuzzificationMethod','centroid');

% ===== Inputs =====
% Temperature ()
fis = addInput(fis,[10 35],'Name','Temperature');
fis = addMF(fis,'Temperature','trapmf',[10 10 15 20],'Name','Cold');      % 1
fis = addMF(fis,'Temperature','trimf' ,[18 22 26]   ,'Name','Comfort');   % 2
fis = addMF(fis,'Temperature','trapmf',[24 28 35 35],'Name','Hot');       % 3

% Light
fis = addInput(fis,[0 100],'Name','LightLevel');
fis = addMF(fis,'LightLevel','trapmf',[0 0 15 35]   ,'Name','Dark');
fis = addMF(fis,'LightLevel','trimf' ,[30 50 70]    ,'Name','Medium');
fis = addMF(fis,'LightLevel','trapmf',[65 85 100 100],'Name','Bright');

% Activity
fis = addInput(fis,[0 1],'Name','Activity');
fis = addMF(fis,'Activity','trapmf',[0 0 0.2 0.4]   ,'Name','Inactive');
fis = addMF(fis,'Activity','trimf' ,[0.3 0.5 0.7]   ,'Name','Moderate');
fis = addMF(fis,'Activity','trapmf',[0.6 0.8 1 1]   ,'Name','Active');

% Humidity
fis = addInput(fis,[0 100],'Name','Humidity');
fis = addMF(fis,'Humidity','trapmf',[0 0 25 35]     ,'Name','Dry');
fis = addMF(fis,'Humidity','trimf' ,[30 50 60]      ,'Name','Comfort');
fis = addMF(fis,'Humidity','trapmf',[55 70 100 100] ,'Name','Humid');

% CO2
fis = addInput(fis,[400 2000],'Name','CO2');
fis = addMF(fis,'CO2','trapmf',[400 400 700 900]    ,'Name','Low');
fis = addMF(fis,'CO2','trimf' ,[800 1100 1400]      ,'Name','Medium');
fis = addMF(fis,'CO2','trapmf',[1300 1600 2000 2000],'Name','High');

% ===== Outputs =====
fis = addOutput(fis,[0 100],'Name','HeaterPower');
fis = addMF(fis,'HeaterPower','trimf',[0 0 30]      ,'Name','Low');    % 1
fis = addMF(fis,'HeaterPower','trimf',[20 40 60]    ,'Name','Medium'); % 2
fis = addMF(fis,'HeaterPower','trimf',[50 100 100]  ,'Name','High');   % 3

fis = addOutput(fis,[0 100],'Name','LightBrightness');
fis = addMF(fis,'LightBrightness','trimf',[0 0 30]  ,'Name','Dim');      % 1
fis = addMF(fis,'LightBrightness','trimf',[25 50 75],'Name','Normal');   % 2
fis = addMF(fis,'LightBrightness','trimf',[70 100 100],'Name','Bright'); % 3

fis = addOutput(fis,[0 100],'Name','BlindsPosition');
fis = addMF(fis,'BlindsPosition','trimf',[0 0 30]   ,'Name','Closed');  % 1
fis = addMF(fis,'BlindsPosition','trimf',[25 50 75] ,'Name','Half');    % 2
fis = addMF(fis,'BlindsPosition','trimf',[70 100 100],'Name','Open');   % 3

fis = addOutput(fis,[0 100],'Name','VentilationFan');
fis = addMF(fis,'VentilationFan','trimf',[0 0 30]   ,'Name','Low');     % 1
fis = addMF(fis,'VentilationFan','trimf',[25 50 75] ,'Name','Medium');  % 2
fis = addMF(fis,'VentilationFan','trimf',[70 100 100],'Name','High');   % 3

%% ===== Rule base =====
% Rule format:
% [Temp Light Act Humid CO2  ->  Heat Light Blind Vent   weight  AND(1)/OR(2)]
comfort_rules = [
  1 0 3 0 0   3 2 2 0   1 1
  1 0 1 0 0   2 1 1 0   1 1
  2 2 2 0 0   1 2 2 0   1 1
  2 3 3 0 0   1 1 1 0   1 1
  3 2 3 0 0   1 1 1 0   1 1
  3 0 1 0 0   1 1 1 0   1 1
  2 2 1 0 0   1 1 1 0   1 1
  1 2 2 0 0   2 2 2 0   1 1
  2 2 3 0 0   1 3 3 0   1 1
  0 2 1 0 0   1 1 1 0   1 1
  0 1 3 0 0   1 3 3 0   1 1
  0 3 1 0 0   1 1 1 0   1 1
];

vent_rules = [
  0 0 0 3 0   0 0 0 3   1 1   % Humid -> Vent High
  0 0 0 0 3   0 0 0 3   1 1   % CO2 High -> Vent High
  0 0 3 0 2   0 0 0 2   1 1   % Active & CO2 Medium -> Vent Medium
  0 0 1 1 1   0 0 0 1   1 1   % Inactive & Dry & Low CO2 -> Vent Low
  0 0 0 2 2   0 0 0 2   1 1   % Comfort humidity & CO2 Medium -> Vent Medium
];

% Important: MATLAB requires at least one nonzero antecedent per rule.
% Therefore, DO NOT add an all-zero "default" antecedent rule here.
fis = addRule(fis, [comfort_rules; vent_rules]);

%% ---------- Save FIS & quick sanity test ----------
writefis(fis, erase(saveAs,'.fis'));
fprintf('Saved FIS -> %s\n', saveAs);

U = [18 20 0.8 60 1400;
     30 80 0.1 40  600;
     22 40 0.5 55  900;
     21 15 0.4 45  850];
Y = evalfis(fis,U);
disp('Sample I/O:');
disp(array2table([U Y], ...
 'VariableNames',{'Temp','Light','Activity','Humidity','CO2','Heater','Brightness','Blinds','Vent'}));

%% ---------- Show rules, MFs, and control surfaces ----------
disp('--- RULES ---');
try
    showrule(fis); % human-readable printout
catch
    R = fis.Rules; for i=1:numel(R), disp(R(i)); end
end

% Rule Viewer
try
    figure('Name','Rule Viewer'); ruleview(fis);
catch
    warning('Rule Viewer not available in this MATLAB session (headless).');
end

% Membership functions
figure('Name','MFs: Inputs'); tiledlayout(2,3);
nexttile; plotmf(fis,'input',1); title('Temperature'); grid on
nexttile; plotmf(fis,'input',2); title('LightLevel'); grid on
nexttile; plotmf(fis,'input',3); title('Activity'); grid on
nexttile; plotmf(fis,'input',4); title('Humidity'); grid on
nexttile; plotmf(fis,'input',5); title('CO2'); grid on

figure('Name','MFs: Outputs'); tiledlayout(2,2);
nexttile; plotmf(fis,'output',1); title('HeaterPower'); grid on
nexttile; plotmf(fis,'output',2); title('LightBrightness'); grid on
nexttile; plotmf(fis,'output',3); title('BlindsPosition'); grid on
nexttile; plotmf(fis,'output',4); title('VentilationFan'); grid on

% Control surfaces (meaningful input pairs)
figure('Name','Control Surfaces'); tiledlayout(2,2);
nexttile; gensurf(fis,[1 3],1); title('Temp & Activity → HeaterPower'); grid on
nexttile; gensurf(fis,[2 3],2); title('Light & Activity → LightBrightness'); grid on
nexttile; gensurf(fis,[2 3],3); title('Light & Activity → BlindsPosition'); grid on
nexttile; gensurf(fis,[5 4],4); title('CO2 & Humidity → VentilationFan'); grid on

% Optional: save PNGs of all open figures
if savePngs
    outdir = fullfile(pwd,'flc_figs'); if ~exist(outdir,'dir'), mkdir(outdir); end
    figs = findall(0,'Type','figure');
    for k=1:numel(figs)
        fname = sprintf('figure_%02d.png', k);
        try
            exportgraphics(figs(k), fullfile(outdir,fname), 'Resolution', 200);
        catch
            saveas(figs(k), fullfile(outdir,fname));
        end
    end
    fprintf('Saved figures to: %s\n', outdir);
end
end
