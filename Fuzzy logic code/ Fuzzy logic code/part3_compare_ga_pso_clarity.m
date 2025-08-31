function part3_compare_ga_pso_clarity()
% PART 3 – GA vs PSO on two CEC'2005 functions: F6 (Shifted Rosenbrock) & F9 (Shifted Rastrigin)
% D = 2 and 10, 15 runs each. Produces clearly labeled CSVs and plots.

clc; close all;

outdir = fullfile(pwd,'part3_outputs_cec');
if ~exist(outdir,'dir'), mkdir(outdir); end
rng(42);  % reproducible seed (state your seed in the report)

% -------------------- Benchmarks (CEC'05) --------------------
benchmarks = {
    struct('id','F6','name','Shifted Rosenbrock', 'lb',-100,'ub',100,'shiftfile','rosenbrock_func_data.mat','bias', 390)
   ,struct('id','F9','name','Shifted Rastrigin',  'lb',  -5,'ub',  5,'shiftfile','rastrigin_func_data.mat', 'bias',-330)
};

dims     = [2 10];
nRuns    = 15;
maxIters = 500;
popSize  = 40;   % GA population / PSO swarm size (same for fairness)

% GA parameters (aligned to Part 2 style)
ga.pc = 0.8; ga.pm = 0.2; ga.tour = 3; ga.mutSigmaFrac = 0.1;

% PSO parameters (classic)
pso.w=0.72; pso.c1=1.49; pso.c2=1.49; pso.vClampFrac=0.2;

% Storage for the combined CSV
rows = {};

% Readme intro
readme = {};
readme{end+1} = 'Part 3 Outputs (GA vs PSO on CEC''05 F6 & F9)';
readme{end+1} = 'Files:';
readme{end+1} = '  - results_part3_cec_summary.csv : Combined summary for all cases';
readme{end+1} = '  - results_F*_D*_GA.csv / results_F*_D*_PSO.csv : Per-case stats (one file per function/dimension/algorithm)';
readme{end+1} = '  - convergence_<FunctionID>_<FunctionName>_D<dim>.png : Convergence (median + IQR) both algorithms';
readme{end+1} = '  - surface_paths_<FunctionID>_<FunctionName>_D2.png   : D=2 surface + contour + median best-paths (GA blue / PSO red)';
readme{end+1} = sprintf('Parameters: Runs=%d, MaxIters=%d, Pop/Swarm=%d, Seed=42', nRuns, maxIters, popSize);
readme{end+1} = ' ';

for b = 1:numel(benchmarks)
    B = benchmarks{b};
    for d = 1:numel(dims)
        D  = dims(d);
        lb = repmat(B.lb,1,D);
        ub = repmat(B.ub,1,D);

        % Load shift vector (if provided), else zeros (still valid)
        o = zeros(1,D);
        if exist(B.shiftfile,'file')
            S = load(B.shiftfile);
            if isfield(S,'o'), o = S.o(1:D); end
        end

        % Function handle wrapper
        if strcmpi(B.id,'F6')
            fwrap = @(x) f6_shifted_rosenbrock(x,o,B.bias);
        else
            fwrap = @(x) f9_shifted_rastrigin(x,o,B.bias);
        end

        % Storage per algorithm for stats/curves
        stats = struct('GA',nan(nRuns,1),'PSO',nan(nRuns,1));
        conv  = struct('GA',nan(nRuns,maxIters),'PSO',nan(nRuns,maxIters));

        % Optional path recording for D=2
        if D==2
            pathGA = nan(nRuns,maxIters,2);
            pathPS = nan(nRuns,maxIters,2);
        end

        % ---------------- GA runs ----------------
        for r=1:nRuns
            if D==2
                [bestFit, bestCurve, bestPath] = run_ga(fwrap,D,lb,ub,popSize,maxIters,ga,true);
                pathGA(r,:,:) = bestPath;
            else
                [bestFit, bestCurve] = run_ga(fwrap,D,lb,ub,popSize,maxIters,ga,false);
            end
            stats.GA(r)   = bestFit;
            conv. GA(r,:) = bestCurve(:)';
            fprintf('[GA ] %s-%s D=%d run %02d | best=%.6e\n',B.id,B.name,D,r,bestFit);
        end

        % ---------------- PSO runs ----------------
        for r=1:nRuns
            if D==2
                [bestFit, bestCurve, bestPath] = run_pso(fwrap,D,lb,ub,popSize,maxIters,pso,true);
                pathPS(r,:,:) = bestPath;
            else
                [bestFit, bestCurve] = run_pso(fwrap,D,lb,ub,popSize,maxIters,pso,false);
            end
            stats.PSO(r)   = bestFit;
            conv. PSO(r,:) = bestCurve(:)';
            fprintf('[PSO] %s-%s D=%d run %02d | best=%.6e\n',B.id,B.name,D,r,bestFit);
        end

        % ---------- Write per-case CSVs (clear names) ----------
        T_GA  = make_stats_table(B.id,B.name,D,'GA', nRuns,maxIters,popSize,stats.GA);
        T_PSO = make_stats_table(B.id,B.name,D,'PSO',nRuns,maxIters,popSize,stats.PSO);

        csv_ga  = fullfile(outdir, sprintf('results_%s_D%d_GA.csv',  B.id, D));
        csv_pso = fullfile(outdir, sprintf('results_%s_D%d_PSO.csv', B.id, D));
        writetable(T_GA,  csv_ga);
        writetable(T_PSO, csv_pso);

        % Add to combined rows
        rows = [rows; table2cell(T_GA); table2cell(T_PSO)];

        % ---------- Convergence: median + IQR shading ----------
        medGA = col_median(conv.GA); medPS = col_median(conv.PSO);
        q1GA  = col_percentile(conv.GA,25); q3GA = col_percentile(conv.GA,75);
        q1PS  = col_percentile(conv.PSO,25); q3PS = col_percentile(conv.PSO,75);

        fig = figure('Color','w','Position',[120 120 900 540]);
        hold on; grid on; set(gca,'YScale','log');
        fill_between(1:maxIters, q1GA, q3GA, [0.75 0.85 1.00], 0.35);
        fill_between(1:maxIters, q1PS, q3PS, [1.00 0.85 0.75], 0.35);
        plot(medGA,'-','LineWidth',1.8,'Color',[0.00 0.20 0.80]);
        plot(medPS,'-','LineWidth',1.8,'Color',[0.80 0.20 0.00]);
        legend({'GA IQR','PSO IQR','GA median','PSO median'},'Location','northeast');
        xlabel('Iteration'); ylabel('Best fitness (log scale)');
        title(sprintf('Convergence – %s (%s)  D=%d', B.id, B.name, D),'FontWeight','bold');
        saveas(fig, fullfile(outdir, sprintf('convergence_%s_%s_D%d.png', B.id, sanitize(B.name), D)));
        close(fig);

        % ---------- For D=2: surface + contour + median paths ----------
        if D==2
            [XX,YY,ZZ,zlo,zhi] = grid_eval_2d(fwrap, lb, ub, 220, 220);
            medPathGA = squeeze(nanmedian(pathGA,1));  % [iter,2]
            medPathPS = squeeze(nanmedian(pathPS,1));  % [iter,2]

            % Surface + paths
            fig = figure('Color','w','Position',[90 90 980 460]);
            tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

            nexttile;
            surf(XX,YY,ZZ,'EdgeColor','none'); view(135,30); grid on; box on;
            colormap(gca,'parula'); caxis([zlo zhi]); colorbar;
            hold on;
            zga = arrayfun(@(i) fwrap(medPathGA(i,:)), 1:size(medPathGA,1))';
            zps = arrayfun(@(i) fwrap(medPathPS(i,:)), 1:size(medPathPS,1))';
            plot3(medPathGA(:,1),medPathGA(:,2),zga,'-','LineWidth',1.9,'Color',[0.00 0.20 0.80]);
            plot3(medPathPS(:,1),medPathPS(:,2),zps,'-','LineWidth',1.9,'Color',[0.80 0.20 0.00]);
            title(sprintf('%s (%s) D=2 Surface + Paths', B.id, B.name),'FontWeight','bold');
            xlabel('x1'); ylabel('x2'); zlabel('f(x)');

            % Contour + paths
            nexttile;
            contourf(XX,YY,ZZ,30,'LineColor','none'); grid on; box on;
            colormap(gca,'parula'); caxis([zlo zhi]); colorbar;
            hold on;
            plot(medPathGA(:,1),medPathGA(:,2),'-','LineWidth',2.0,'Color',[0.00 0.20 0.80]);
            plot(medPathPS(:,1),medPathPS(:,2),'-','LineWidth',2.0,'Color',[0.80 0.20 0.00]);
            legend({'','GA path','PSO path'},'Location','northeast');
            title(sprintf('%s (%s) D=2 Contours + Paths', B.id, B.name),'FontWeight','bold');
            xlabel('x1'); ylabel('x2');

            saveas(fig, fullfile(outdir, sprintf('surface_paths_%s_%s_D2.png', B.id, sanitize(B.name))));
            close(fig);
        end

        % Update README lines
        readme{end+1} = sprintf('%s (%s), D=%d:', B.id, B.name, D);
        readme{end+1} = sprintf('  - %s', rel(outdir, sprintf('results_%s_D%d_GA.csv',  B.id, D)));
        readme{end+1} = sprintf('  - %s', rel(outdir, sprintf('results_%s_D%d_PSO.csv', B.id, D)));
        readme{end+1} = sprintf('  - %s', rel(outdir, sprintf('convergence_%s_%s_D%d.png', B.id, sanitize(B.name), D)));
        if D==2
            readme{end+1} = sprintf('  - %s', rel(outdir, sprintf('surface_paths_%s_%s_D2.png', B.id, sanitize(B.name))));
        end
        readme{end+1} = ' ';
    end
end

% ---------- Combined CSV ----------
Summary = cell2table(rows, 'VariableNames', ...
    {'FunctionID','FunctionName','Dim','Algorithm','Runs','MaxIters','PopOrSwarm','Mean','Std','Best','Worst'});
writetable(Summary, fullfile(outdir,'results_part3_cec_summary.csv'));

% ---------- MAT dump ----------
save(fullfile(outdir,'results_part3_cec.mat'), 'Summary','benchmarks','dims','nRuns','maxIters','popSize','ga','pso');

% ---------- README ----------
fid = fopen(fullfile(outdir,'README.txt'),'w');
for i=1:numel(readme), fprintf(fid,'%s\n',readme{i}); end
fclose(fid);

% ---------- Pretty console print ----------
disp('=================== PART 3 SUMMARY (clearly labeled) ===================');
disp(Summary);
disp(['Saved outputs in: ', outdir]);

end % main

% ========================= Helper: per-case table =========================
function T = make_stats_table(fid,fname,D,algo,nRuns,maxIters,pop,vals)
    T = table( string(fid), string(fname), D, string(algo), nRuns, maxIters, pop, ...
               mean(vals), std(vals), min(vals), max(vals), ...
               'VariableNames', {'FunctionID','FunctionName','Dim','Algorithm','Runs','MaxIters','PopOrSwarm','Mean','Std','Best','Worst'});
end

% ============================= Math helpers ===============================
function m = col_median(M),  m = prctile(M,50,1); end
function q = col_percentile(M,p), q = prctile(M,p,1); end

% ============================ Plotting helper =============================
function fill_between(x,lo,hi,RGB,alphaVal)
    xx = [x(:).' fliplr(x(:).')];
    yy = [lo(:).' fliplr(hi(:).')];
    p  = patch(xx,yy,RGB,'EdgeColor','none'); set(p,'FaceAlpha',alphaVal);
end

% ================================ Surfaces ================================
function [XX,YY,ZZ,zlo,zhi] = grid_eval_2d(fwrap, lb, ub, nx, ny)
    xs = linspace(lb(1),ub(1),nx); ys = linspace(lb(2),ub(2),ny);
    [XX,YY] = meshgrid(xs,ys);
    ZZ = zeros(size(XX));
    for ii = 1:numel(XX)
        ZZ(ii) = fwrap([XX(ii) YY(ii)]);
    end
    zlo = prctile(ZZ(:),1); zhi = prctile(ZZ(:),99);
end

% ============================ GA / PSO cores ==============================
function [bestFit,bestCurve,bestPath] = run_ga(f,D,lb,ub,popSize,maxIters,ga,wantPath)
    if nargin<8, wantPath=false; end
    if wantPath && D==2, bestPath = nan(maxIters,2); else, bestPath = []; end
    P = rand(popSize,D).*(ub-lb)+lb;
    fit = arrayfun(@(i) f(P(i,:)), 1:popSize)'; bestCurve = nan(maxIters,1);
    bestFit = inf; bestX = zeros(1,D);
    for it=1:maxIters
        [fbest, ibest] = min(fit);
        if fbest < bestFit, bestFit=fbest; bestX=P(ibest,:); end
        bestCurve(it) = bestFit;
        if ~isempty(bestPath), bestPath(it,:) = bestX; end
        Q = zeros(size(P));
        for k=1:2:popSize
            p1=tournament(fit,ga.tour); p2=tournament(fit,ga.tour);
            c1=P(p1,:); c2=P(p2,:);
            if rand<ga.pc
                cx=randi(D-1); tmp=c1(cx+1:end); c1(cx+1:end)=c2(cx+1:end); c2(cx+1:end)=tmp;
            end
            sig=ga.mutSigmaFrac*(ub-lb);
            if rand<ga.pm, c1=c1+sig.*randn(1,D); end
            if rand<ga.pm, c2=c2+sig.*randn(1,D); end
            c1=min(ub,max(lb,c1)); c2=min(ub,max(lb,c2));
            Q(k,:)=c1; if k+1<=popSize, Q(k+1,:)=c2; end
        end
        [~,iworst]=max(fit); Q(iworst,:)=bestX;
        P=Q; fit=arrayfun(@(i) f(P(i,:)), 1:popSize)';
    end
end
function idx=tournament(fit,t), cand=randi(numel(fit),[t,1]); [~,k]=min(fit(cand)); idx=cand(k); end

function [bestFit,bestCurve,bestPath] = run_pso(f,D,lb,ub,swarmSize,maxIters,p,wantPath)
    if nargin<8, wantPath=false; end
    if wantPath && D==2, bestPath = nan(maxIters,2); else, bestPath = []; end
    X=rand(swarmSize,D).*(ub-lb)+lb; V=zeros(swarmSize,D);
    pBest=X; pBestFit=arrayfun(@(i) f(X(i,:)), 1:swarmSize)'; [gBestFit,gIdx]=min(pBestFit); gBest=pBest(gIdx,:);
    bestCurve=nan(maxIters,1); vClamp=p.vClampFrac*(ub-lb);
    for it=1:maxIters
        r1=rand(swarmSize,D); r2=rand(swarmSize,D);
        V=p.w*V+p.c1*r1.*(pBest-X)+p.c2*r2.*(repmat(gBest,swarmSize,1)-X);
        V=sign(V).*min(abs(V),repmat(vClamp,swarmSize,1));
        X=X+V; X=min(ub,max(lb,X));
        fit=arrayfun(@(i) f(X(i,:)), 1:swarmSize)';
        improved=fit<pBestFit; pBest(improved,:)=X(improved,:); pBestFit(improved)=fit(improved);
        [currBest,idx]=min(pBestFit); if currBest<gBestFit, gBestFit=currBest; gBest=pBest(idx,:); end
        bestCurve(it)=gBestFit;
        if ~isempty(bestPath), bestPath(it,:)=gBest; end
    end
    bestFit=gBestFit;
end

% ============================ CEC'05 functions ============================
function val = f6_shifted_rosenbrock(x,o,bias)
    z = x - o(:)'; D = numel(z);
    if D<2, val = (z(1)-1)^2 + bias; return; end
    s=0;
    for i=1:(D-1)
        s = s + 100*(z(i+1)-z(i)^2)^2 + (z(i)-1)^2;
    end
    val = s + bias;
end
function val = f9_shifted_rastrigin(x,o,bias)
    z = x - o(:)'; val = sum(z.^2 - 10*cos(2*pi*z) + 10) + bias;
end

% ============================ Utils ============================
function s = sanitize(str), s = regexprep(char(str),'[^A-Za-z0-9_-]','_'); end
function r = rel(folder, fname), r = fullfile(folder, fname); end

