function fig3A(analDir,rrmDir,tsDir,figDir)
% function Fig3A(analDir,rrmDir,tsDir,figDir)
%
% plot Fig3A: the spike-density-functions averaged across all the sessions, 
% separated by movement and attention condition
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

% data SDFs
bsdfsStM = cell(3,2);
bsdfsStNM = cell(3,2);

% cross-validated full model
bsdfsVfStM = cell(3,2);
bsdfsVfStNM = cell(3,2);

for n = 1:length(dn)
    sessionID = dn(n).name;
    
    % load data
    load(fullfile(rrmDir,sessionID,'analysisVars.mat'),'Vc');
    load(fullfile(tsDir,sessionID,'preprocessed_data'),'D');
    load(fullfile(analDir,sessionID,'sdfData'),'sdfStM','sdfStNMM','sdfMM','sdfNMM');
    load(fullfile(analDir,sessionID,'sdfFullModel'),'sdfFModStM','sdfFModStNMM','sdfFModMM','sdfFModNMM','tSDF');    
    
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

                % sdfs for averaging                
                bsdfsStM{v,1} = [bsdfsStM{v,1};cell2mat(sdfStM(idx,1))];
                bsdfsStM{v,2} = [bsdfsStM{v,2};cell2mat(sdfStM(idx,2))];
                bsdfsStNM{v,1} = [bsdfsStNM{v,1};cell2mat(sdfStNMM(idx,1))];
                bsdfsStNM{v,2} = [bsdfsStNM{v,2};cell2mat(sdfStNMM(idx,2))];                

                bsdfsVfStM{v,1} = [bsdfsVfStM{v,1};cell2mat(sdfFModStM(idx,1))];
                bsdfsVfStM{v,2} = [bsdfsVfStM{v,2};cell2mat(sdfFModStM(idx,2))];
                bsdfsVfStNM{v,1} = [bsdfsVfStNM{v,1};cell2mat(sdfFModStNMM(idx,1))];
                bsdfsVfStNM{v,2} = [bsdfsVfStNM{v,2};cell2mat(sdfFModStNMM(idx,2))];
        end        
    end % area
end % session


%%
fh=figure;
pos=get(fh,'position');
set(fh,'position',[pos(1:2) 415 305])
icolM =  flipud([.6 .6 1; 0 0 1]);
icolNM = flipud([1 .8 .3; 1 .5 .1 ]);

kwidth = 2; % SD of smoothing kernel width in ms
x = 0:kwidth*3;
half_g = exp(-(x.^2)/(2 * kwidth^2));
scale = 60/sum(half_g); % ms resolution for spikes
dt  = length(half_g)/2;  %

% data sdfs
for n = 1
    sp{6}=subplot(1,2,1);
    
        hold on
    g=plot(tSDF+dt,conv(nanmean([...
        cell2mat(bsdfsStM(1,2));...
        cell2mat(bsdfsStM(2,2));...
        cell2mat(bsdfsStM(3,2));...
        ]),half_g,'same')*scale,'-b','linewidth',2);
    g.Color = icolM(2,:);

    
    g=plot(tSDF+dt,conv(nanmean([...
        cell2mat(bsdfsStM(1,1));...
        cell2mat(bsdfsStM(2,1));...
        cell2mat(bsdfsStM(3,1))...
        ]),half_g,'same')*scale,'-b','linewidth',2);
    g.Color = icolM(1,:);
    
            g=plot(tSDF+dt,conv(nanmean([...
        cell2mat(bsdfsStNM(1,2));...
        cell2mat(bsdfsStNM(2,2));...
        cell2mat(bsdfsStNM(3,2))...
        ]),half_g,'same')*scale,'-b','linewidth',1.5);
    g.Color = icolNM(2,:);


    
    g=plot(tSDF+dt,conv(nanmean([...
        cell2mat(bsdfsStNM(1,1));...
        cell2mat(bsdfsStNM(2,1));...
        cell2mat(bsdfsStNM(3,1))...
        ]),half_g,'same')*scale,'-b','linewidth',1.5);
    g.Color = icolNM(1,:);
    
    
    if n ==1
        xlabel('time after stimulus')
        ylabel('response [sp/sec]');
    end
offsetAxes;
end

text(1250, 57,'attention')
text(1290, 52,'in   out')
text(200, 50,'with motion')
text(200, 45,'w/o motion')
plot([1300 1500],ones(1,2)*50,'-','linewidth',2,'color',icolM(1,:));
plot([1600 1800],ones(1,2)*50,'-','linewidth',2,'color',icolM(2,:));
plot([1300 1500],ones(1,2)*45,'-','linewidth',2,'color',icolNM(1,:));
plot([1600 1800],ones(1,2)*45,'-','linewidth',2,'color',icolNM(2,:));


% model SDFs
for n = 4
    sp{7}=subplot(1,2,2);
    hold on
    g=plot(tSDF+dt,conv(nanmean([...
        cell2mat(bsdfsVfStM(1,2));...
        cell2mat(bsdfsVfStM(2,2));...
        cell2mat(bsdfsVfStM(3,2))...
        ]),half_g,'same')*scale,'-b','linewidth',2);
    g.Color = icolM(2,:);
        g=plot(tSDF+dt,conv(nanmean([...
        cell2mat(bsdfsVfStM(1,1));...
        cell2mat(bsdfsVfStM(2,1));...
        cell2mat(bsdfsVfStM(3,1))...
        ]),half_g,'same')*scale,'-b','linewidth',2);
    g.Color = icolM(1,:);

    
    
        g=plot(tSDF+dt,conv(nanmean([...
        cell2mat(bsdfsVfStNM(1,2));...
        cell2mat(bsdfsVfStNM(2,2));...
        cell2mat(bsdfsVfStNM(3,2))...
        ]),half_g,'same')*scale,'-b','linewidth',1.5);
    g.Color = icolNM(2,:);
    
        g=plot(tSDF+dt,conv(nanmean([...
        cell2mat(bsdfsVfStNM(1,1));...
        cell2mat(bsdfsVfStNM(2,1));...
        cell2mat(bsdfsVfStNM(3,1))...
        ]),half_g,'same')*scale,'-b','linewidth',1);
    g.Color = icolNM(1,:);

    
    if n ==1
        xlabel('time after stimulus (ms)')
        %ylabel('response [sp/sec]');
    end
    offsetAxes;
    
end
for n = 6:7
    set(sp{n},'xlim',[-200 1900],'xtick', 0:600:1800, 'box','off','tickdir','out','ylim',[0 60],'ytick',[0 30 60]);
end

if ~exist(figDir,'dir')
    mkdir(figDir)
end
cd(figDir)
fname = ['fig3A_SDFs' date '.pdf'];
figname = ['fig3A_SDFs' date '.fig'];
saveas(gcf,fname)
saveas(gcf,figname)

% save figure data 
datname = ['fig3A_SDF_data' date '.mat'];
save(datname,'bsdfsStM','bsdfsStNM', 'bsdfsVfStM', 'bsdfsVfStNM');
