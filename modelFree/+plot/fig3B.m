function fig3B(analDir,rrmDir,tsDir,figDir)
% function Fig3b(analDir,rrmDir,tsDir,figDir)
%
% plot Fig3b: Movement Index and Attention Index for each unit
% save figure as *fig and pdf
%
% input:
% analDir   path to directory of results for model-free analysis
% rrmDir    path to directory of preprocessed data for RRM fitting
%           to compute inclusion criterion and get informrmation on visual
%           area
% tsDir     path to directory of results from RRM fitting to compute
%           inclusion criterion
% figDir    directory in which to save the figures

% get list of sessions for M2
cd(analDir)
dn = dir('M2_*'); 

%% harvest data for plotting ==============================================
wMIs_Vdt = cell(3,1);
wAIs_Vdt = cell(3,1);

for n = 1:length(dn)

    sessionID = dn(n).name;
    
    % load data
    load(fullfile(rrmDir,sessionID,'analysisVars.mat'),'Vc');
    load(fullfile(tsDir,sessionID,'preprocessed_data'),'D');
    load(fullfile(analDir,sessionID,'AI_MIs'),'AI_Vdt','wMI_Vdt')
    
    % inclusion criterion based on mean rates
    [~,mr ]         = io.stim_rates_ma(Vc,D);
    mrb             = sum(mr>2,2)';

    for v = 1:3
        area_idx = [];
        if isfield(D.area,['V' num2str(v) '_idx'])
            area_idx = eval(['D.area.V' num2str(v) '_idx;']);
        end
        
        % select unit in respective area exceeding required rates
        idx = intersect(area_idx, find(mrb==4));
        if ~isempty(idx) 

            wMIs_Vdt{v} = [wMIs_Vdt{v}, wMI_Vdt(idx)];
            wAIs_Vdt{v} = [wAIs_Vdt{v}, AI_Vdt(idx)];
        end
    end
end

%% plot data ==============================================================
fh=figure;
pos = get(fh,'position');
set(fh,'position',[pos(1:2)/2 370 350])

% scatterplot
wi = .4;
hi = 0.4;
bo = .15;
le = 0.15;

set(fh,'name','AIvsMI by area')
subplot('position',[le bo wi hi])
plot([-.35 .35],zeros(1,2),'--','linewidth',1,'color',[.5 .5 .5]);
hold on;
plot(zeros(1,2),[-.35 .35],'--','linewidth',1,'color',[.5 .5 .5]);
plot([-.35 .35],[-.35 .35],'-','linewidth',1,'color',[.5 .5 .5])

xlims = [-.63 :0.05:.63];
vcols = [0.3804    0.2510    0.5804; ...
    0.0549    0.4588    0.5608; ...
    0.7490    0.8118    0.5490];

for n = 1:3
    g=plot(wAIs_Vdt{4-n},wMIs_Vdt{4-n},'o');
    g.MarkerFaceColor = vcols(4-n,:);
    g.Color = vcols(4-n,:);
    g.MarkerSize = 3;
end
set(gca,'xlim',[-.35 .35],'ylim',[-.35 .35],'box','off','tickdir','out',...
    'xtick',[-.3 0 .3],'ytick',[-.3 0 .3]);
offsetAxes;
xlabel('attention index')
ylabel('movement index')


% histograms for AI
wi = .4;
hi = 0.08;
xhi = .1;
bo = .6;
le = 0.15;
xlims = [-.65 :0.05:.65];
for n = 1:3
    subplot('position',[le, bo+(n-1)*xhi, wi hi])
    
    hist(wAIs_Vdt{n},xlims);
    
    
    oh = get(gca,'children');
    oh(1).FaceColor=vcols(n,:);
    oh(1).EdgeColor = 'none';
    
    hold on;
    yl = get(gca,'ylim');
    plot(zeros(1,2),yl*1.2,'--','linewidth',1,'color',[.5 .5 .5]);
    
    h=plot(mean(wAIs_Vdt{n}),yl(2)*1.15,'v');
    h.Color=vcols(n,:);
    h.MarkerEdgeColor = vcols(n,:);
    h.MarkerFaceColor = vcols(n,:);
    set(gca,'xlim',[-.35 .35],'box','off','tickdir','out',...
        'xtick',[-.3 0 .3],'xticklabel','', 'ytick',yl,'yticklabel',yl);
    text(-.3,yl(2),sprintf('n=%d',sum(~isnan(wAIs_Vdt{n}))));
    offsetAxes;
end

% histograms for MI
hi = .4;
wi = 0.08;
xwi = .1;
bo = .15;
le = 0.6;
for n = 1:3
    
    subplot('position',[le+ (n-1)*xwi, bo wi hi])
    
    hist(wMIs_Vdt{n},xlims);
    oh = get(gca,'children');
    oh(1).FaceColor=vcols(n,:);
    oh(1).EdgeColor = 'none';
    
    hold on;
    yl = get(gca,'ylim');
    plot(zeros(1,2),yl*1.2,'--','linewidth',1,'color',[.5 .5 .5]);
    
    h=plot(mean(wMIs_Vdt{n}),yl(2)*1.15,'<');
    h.Color=vcols(n,:);
    h.MarkerEdgeColor = vcols(n,:);
    h.MarkerFaceColor = vcols(n,:);
    
    set(gca,'tickdir','out','xtick',[-.3 0 .3],'xticklabel','',...
         'ytick',yl,'yticklabel',yl);
    set(gca,'xlim',[-.35 .35],'box','off','xdir','normal','view',[90 -90]);
    offsetAxes;
end

% save figure =============================================================
if ~exist(figDir,'dir')
    mkdir(figDir)
end
cd(figDir)

fname = ['fig3B_AIvsMI' date '.pdf'];
figname = ['fig3B_AIvsMI' date '.fig'];
saveas(gcf,fname)
saveas(gcf,figname)

% save figure data 
datname = ['fig3B_AIvsMI_data' date '.mat'];
save(datname,'wMIs_Vdt','wAIs_Vdt')
