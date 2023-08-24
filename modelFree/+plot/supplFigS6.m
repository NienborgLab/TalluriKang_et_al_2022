function supplFigS6(analDir,rrmDir,tsDir,figDir)
% function supplFig6(analDir,rrmDir,tsDir,figDir)
%
% plot supplFig6: PSTHs across sessions separately for M1 and M2, Movement
% index across all units.
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
%


%%
cd(analDir)
dn = dir('*_*');

% initialize
wMIs_Vdt = cell(3,1);

%data SDFs
bsdfsStM = cell(3,2);
bsdfsStNM = cell(3,2);
bsdfsM = cell(3,2);
bsdfsNM = cell(3,2);


tSDFs = cell(1,2);


for n = 1:length(dn)
        sessionID = dn(n).name;

    % load data
    load(fullfile(rrmDir,sessionID,'analysisVars.mat'),'Vc');
    load(fullfile(tsDir,sessionID,'preprocessed_data'),'D');

    if contains(sessionID,'M1')
        [~,mr ]         = io.stim_rates_le(Vc,D);
        load(fullfile(analDir,sessionID,'MIs'),'wMI_Vdt');
        load(fullfile(analDir,sessionID,'sdfData'),'sdfStM',...
            'sdfStNMM','sdfMM','sdfNMM','tSDF');
        tSDFs{1} = tSDF;
    else
        [~,mr ]         = io.stim_rates_ma(Vc,D);
        load(fullfile(analDir,sessionID,'AI_MIs'),'wMI_Vdt')
        load(fullfile(analDir,sessionID,'sdfData'),'sdfStM',...
            'sdfStNMM','sdfMM','sdfNMM','tSDF');
        tSDFs{2} = tSDF;
    end

    % inclusion criterion based on mean rates
    mrb             = sum(mr>2,2)';

    for v = 1:3 % visual areas
        area_idx = [];
        if isfield(D.area,['V' num2str(v) '_idx'])
            area_idx = eval(['D.area.V' num2str(v) '_idx;']);
        end

        % select unit in respective area exceeding required rates
        idx = intersect(area_idx, find(mrb==4));

        if ~isempty(idx) 

            wMIs_Vdt{v} = [wMIs_Vdt{v}, wMI_Vdt(idx)];

            if contains(sessionID,'M1')
                % sdfs for averaging only high contrast stimuli
                bsdfsM{v,1} = [bsdfsM{v,1};cell2mat(sdfMM(idx,4))];
                bsdfsNM{v,1} = [bsdfsNM{v,1};cell2mat(sdfNMM(idx,4))];

                bsdfsStM{v,1} = [bsdfsStM{v,1};cell2mat(sdfStM(idx,4))];
                bsdfsStNM{v,1} = [bsdfsStNM{v,1};cell2mat(sdfStNMM(idx,4))];

            else
                % sdfs for averaging
                bsdfsM{v,2} = [bsdfsM{v,2};cell2mat(sdfMM(idx,1));cell2mat(sdfMM(idx,2))];
                bsdfsNM{v,2} = [bsdfsNM{v,2};cell2mat(sdfNMM(idx,1));cell2mat(sdfNMM(idx,2))];

                bsdfsStM{v,2} = [bsdfsStM{v,2};cell2mat(sdfStM(idx,1));cell2mat(sdfStM(idx,2))];
                bsdfsStNM{v,2} = [bsdfsStNM{v,2};cell2mat(sdfStNMM(idx,1));cell2mat(sdfStNMM(idx,2))];
            end                    
        end            
    end  % visual area  
end % list

%% plot data
figure
set(gcf,'position',[61   286   889   509])
wi = .14;
hi = 0.21;
bo = .72;
le = [0.1 .45];
icolM =  flipud([.6 .6 1; 0 0 1]);
icolNM = flipud([1 .8 .3; 1 .5 .1 ]);
xlimsM = [-450 450];
xlims = [-250 450];
ylims = NaN(1,2);

fS = 60; % sampling frequency (Hz) of data
kwidth = 1; % in 16ms time-bin 
x = 0:kwidth*3;
half_g = exp(-(x.^2)/(2 * kwidth^2));
scale = 1/sum(half_g); % y-rescaling
dT    = 24; % in ms (ensure SDFs are centered on 0)
scl = [fS, fS]*scale; 

 
% sdfs triggered on stimulus onset
for n = 1:2
    subplot('position',[le(n),bo,wi,hi])
    hold on
    g=plot(tSDFs{n}+dT,scl(n)*conv(nanmean([...
        cell2mat(bsdfsStM(1,n));...
        cell2mat(bsdfsStM(2,n));...
        cell2mat(bsdfsStM(3,n));...
        ]),half_g,'same'),'-b','linewidth',2);
    g.Color = icolM(1,:);

            g=plot(tSDFs{n}+dT,scl(n)*conv(nanmean([...
        cell2mat(bsdfsStNM(1,n));...
        cell2mat(bsdfsStNM(2,n));...
        cell2mat(bsdfsStNM(3,n))...
        ]),half_g,'same'),'-b','linewidth',1.5);
    g.Color = icolNM(1,:);

    xlabel('   Time [ms]\newlineafter Stimulus')
    if n==1
        ylabel('Response [sp/sec]')
    end   
    yls = get(gca,'ylim');
    ylims(n) = yls(2);
    set(gca,'xlim',xlims,'ylim',[0 yls(2)],'xtick',[-300:150:450],'tickdir','out')
    title(sprintf...
        ('%s %d','                                             Data: M',n))
    offsetAxes;
    if n ==1
        text(-700,yls(2)*1.2,'A','fontsize',16,'fontweight','bold')
    else
        text(-700,yls(2)*1.2,'B','fontsize',16,'fontweight','bold')
    end
   
end 


% sdfs triggered on motion onset ==========================================
le = [0.26 .61];
for n = 1:2
    
    subplot('position',[le(n),bo, wi,hi])
    hold on
    g=plot(tSDFs{n}+dT,scl(n)*conv(nanmean([...
        cell2mat(bsdfsM(1,n));...
        cell2mat(bsdfsM(2,n));...
        cell2mat(bsdfsM(3,n));...
        ]),half_g,'same'),'-b','linewidth',2);
    g.Color = icolM(1,:);

            g=plot(tSDFs{n}+dT,scl(n)*conv(nanmean([...
        cell2mat(bsdfsNM(1,n));...
        cell2mat(bsdfsNM(2,n));...
        cell2mat(bsdfsNM(3,n))...
        ]),half_g,'same'),'-b','linewidth',1.5);
    g.Color = icolNM(1,:);

    xlabel('\newlineafter Movement')
    yls = get(gca,'ylim');
    set(gca,'xlim',xlimsM,'ylim',[0 ylims(n)],'yticklabel','','xtick',...
        [-300:150:450],'tickdir','out')
    offsetAxes;
    if n==1
        %%
        text(-200, 28,'movement:')
        text(-50, 24,'with')
        text(-50, 20,'w/o')
        plot([200 400],ones(1,2)*24,'-','linewidth',2,'color',icolM(1,:));
        plot([200 400],ones(1,2)*20,'-','linewidth',2,'color',icolNM(1,:));
        
    end
    
end 



% histograms for MI =======================================================
wi = .13;
hi = 0.062;
xhi = .072;
bo = 0.72;
le = 0.8;

xlims = [-.65 :0.05:.65];
vcols = [0.3804    0.2510    0.5804; ...
    0.0549    0.4588    0.5608; ...
    0.7490    0.8118    0.5490];

for n = 1:3
    subplot('position',[le, bo+(n-1)*xhi, wi hi])
    
    hist(wMIs_Vdt{n},xlims);
    oh = get(gca,'children');
    oh(1).FaceColor=vcols(n,:);
    oh(1).EdgeColor = 'none';
    hold on;
    yl = get(gca,'ylim');
    plot(zeros(1,2),yl*1.2,'--','linewidth',1,'color',[.5 .5 .5]);
    
    h=plot(mean(wMIs_Vdt{n}),yl(2)*1.2,'v');
    h.Color=vcols(n,:);
    h.MarkerEdgeColor = vcols(n,:);
    h.MarkerFaceColor = vcols(n,:);
    set(gca,'xlim',[-.35 .35],'box','off','tickdir','out',...
        'xtick',[-.3 0 .3],'xticklabel','', 'ytick',yl,'yticklabel',yl);
    offsetAxes;
    text(-.32,yl(2),sprintf('n=%d',sum(~isnan(wMIs_Vdt{n}))))
    if n ==1
        xlabel('       Movement Index\newline (R_{with}-R_{w/o})/(R_{with}+R_{w/o})')
        set(gca,'xticklabel',[-.3 0 .3],'yticklabel',[0:100:500])
    elseif n==2
        ylabel ('   # Units\newline')
    else
        text(-0.5,400,'C','fontsize',16,'fontweight','bold')
    end
end

if ~exist(figDir,'dir')
    mkdir(figDir)
end

% save figure =============================================================
cd(figDir)
fname = ['MI_SDFs' date '.pdf'];
figname = ['MI_SDFs' date '.fig'];
saveas(gcf,fname)
saveas(gcf,figname)

% save figure data 
datname = ['figS6_MI_data' date '.mat'];
save(datname,'wMIs_Vdt')
