%% Part 2 – GA Optimisation of FLC Membership Functions 

clear; clc; close all;

% --- Load baseline FIS from Part 1
fis0 = readfis('AssistiveHomeFLC_Ext.fis');

% --- Ensure coverage: "Any" MF + soft-default rule
fis0 = p2_ensure_soft_default_coverage(fis0);

% --- dataset 
N = 500;
T = 10 + (35-10)*rand(N,1);
L = 100*rand(N,1);
A = rand(N,1);
H = 100*rand(N,1);
C = 400 + (2000-400)*rand(N,1);
X = [T L A H C];
Y0 = evalfis(fis0,X);
Y  = max(0,min(100, Y0 + randn(size(Y0)).*[3 3 3 3]));

% --- GA Parameters
popSize = 30; maxGen = 60; pc = 0.8; pm = 0.2; tour=3;
bounds  = [10 35; 0 100; 0 1; 0 100; 400 2000; 0 100; 0 100; 0 100; 0 100];
nVars   = size(bounds,1)*6;  % 6 genes per variable

% --- Initialise population
P = rand(popSize,nVars); fit = inf(popSize,1);
best = inf; bestG = []; history = nan(maxGen,1);

for g = 1:maxGen
    for i=1:popSize
        if ~isfinite(fit(i))
            fit(i) = p2_fitness(P(i,:),fis0,bounds,X,Y);
        end
    end
    [fg,ix] = min(fit);
    if fg < best, best=fg; bestG=P(ix,:); end
    history(g)=best;
    fprintf('Gen %d | Best MSE=%.4f\n',g,best);

    % --- Reproduction
    Q = zeros(size(P));
    for k=1:2:popSize
        p1=p2_tournament(fit,tour); p2_=p2_tournament(fit,tour);
        c1=P(p1,:); c2=P(p2_,:);
        if rand<pc
            cx=randi(nVars-1);
            tmp=c1(cx+1:end); c1(cx+1:end)=c2(cx+1:end); c2(cx+1:end)=tmp;
        end
        if rand<pm, c1=c1+0.05*randn(1,nVars); end
        if rand<pm, c2=c2+0.05*randn(1,nVars); end
        Q(k,:)=min(1,max(0,c1));
        if k+1<=popSize, Q(k+1,:)=min(1,max(0,c2)); end
    end
    [~,w]=max(fit); Q(w,:)=bestG;  % elitism
    P=Q; fit(:)=inf;
end

% --- Build and save best FIS
fisBest = p2_rebuild(bestG,fis0,bounds);
writefis(fisBest,'AssistiveHomeFLC_Ext_GA');
fprintf('Optimised FIS saved as AssistiveHomeFLC_Ext_GA.fis\n');

% --- Plot convergence
figure; plot(history,'LineWidth',1.5); grid on;
xlabel('Generation'); ylabel('Best MSE'); title('GA Optimisation Convergence');

%% ===== Helper functions =====
function idx=p2_tournament(f,t)
    cand=randi(numel(f),[t,1]); [~,k]=min(f(cand)); idx=cand(k);
end

function mse=p2_fitness(genes,fis0,bounds,X,Y)
    fis=p2_rebuild(genes,fis0,bounds);
    Yhat=evalfis(fis,X);
    mse=mean((Yhat(:)-Y(:)).^2,'omitnan');
    % mild penalty if a mid-range fallback slips in
    if any(abs(Yhat(:)-50)<1e-9), mse=mse+1; end
end

function fis=p2_rebuild(genes,fis0,bounds)
    % Update MF parameters IN PLACE to keep rule indices intact
    fis=fis0; per=6;
    % INPUT indices: 1=Temp,2=Light,3=Activity,4=Humidity,5=CO2
    % OUTPUT indices: 1=Heater,2=Brightness,3=Blinds,4=Vent
    for v=1:size(bounds,1)
        g=genes((v-1)*per+(1:per));
        k=p2_knots(g,bounds(v,:));
        if v==1
            % Temperature: tune MF #1..#3 (types must match Part 1)
            fis = setMFParamsInPlace(fis,true ,1,1,[k(1) k(1) k(2) k(3)]); % trapmf
            fis = setMFParamsInPlace(fis,true ,1,2,[k(2) k(3) k(4)]);       % trimf
            fis = setMFParamsInPlace(fis,true ,1,3,[k(4) k(5) k(6) k(6)]); % trapmf
            % MF #4 'Any' untouched
        elseif v<=5
            inIdx=v;
            fis = setMFParamsInPlace(fis,true ,inIdx,1,[k(1) k(1) k(2) k(3)]); % trapmf
            fis = setMFParamsInPlace(fis,true ,inIdx,2,[k(2) k(3) k(4)]);       % trimf
            fis = setMFParamsInPlace(fis,true ,inIdx,3,[k(4) k(5) k(6) k(6)]); % trapmf
        else
            outIdx=v-5;
            fis = setMFParamsInPlace(fis,false,outIdx,1,[k(1) k(1) k(3)]);     % trimf
            fis = setMFParamsInPlace(fis,false,outIdx,2,[k(2) k(3) k(4)]);     % trimf
            fis = setMFParamsInPlace(fis,false,outIdx,3,[k(4) k(6) k(6)]);     % trimf
        end
    end
end

function k=p2_knots(u,dom)
    u=max(1e-3,min(1-1e-3,u(:))); w=u/sum(u); c=cumsum(w);
    k=dom(1)+c'*(dom(2)-dom(1));
    for i=2:numel(k), if k(i)<=k(i-1), k(i)=k(i-1)+1e-4*(dom(2)-dom(1)); end; end
end

% ---------- NEW: set params IN PLACE (no remove/add) ----------
function fis=setMFParamsInPlace(fis,isInput,varIdx,mfIdx,params)
    if isInput
        % keep MF type as created in Part 1 (trap/tri)
        fis.Inputs(varIdx).MembershipFunctions(mfIdx).Parameters = params;
    else
        fis.Outputs(varIdx).MembershipFunctions(mfIdx).Parameters = params;
    end
end

function fis=p2_ensure_soft_default_coverage(fis)
    % Ensure Temperature has "Any" MF and default rule to Low outputs
    TR=fis.Inputs(1).Range;
    mfNames=string({fis.Inputs(1).MembershipFunctions.Name});
    if ~any(lower(mfNames)=="any")
        fis=addMF(fis,'Temperature','trapmf',[TR(1) TR(1) TR(2) TR(2)],'Name','Any');
    end
    mfNames=string({fis.Inputs(1).MembershipFunctions.Name});
    anyIdx=find(lower(mfNames)=="any",1); needDefault=true;
    for r=1:numel(fis.Rules)
        ant=fis.Rules(r).Antecedent; con=fis.Rules(r).Consequent;
        if ant(1)==anyIdx && all(ant(2:end)==0) && all(con(1:4)==1), needDefault=false; end
    end
    if needDefault, fis=addRule(fis,[anyIdx 0 0 0 0 1 1 1 1 0.3 1]); end
end

