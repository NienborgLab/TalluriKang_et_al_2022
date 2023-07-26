function plotOneEx(tbl)

thr = 0.8;

nCh = size(tbl,1);

switch nCh
    case 1
        nrow = 1;
        ncol = 1;
    case 24
        nrow = 6+1;
        ncol = 4;
        ha_ch = gobjects(4,6);
    case 32
        nrow = 8+1;
        ncol = 4;
        ha_ch = gobjects(4,8);
    otherwise
end

lcolor = lines(2);
lineWidth = 1.5;

clf
xx = linspace(tbl.fit{1}.dat.x(1),tbl.fit{1}.dat.x(end),100);
xlim = tbl.fit{1}.dat.x([1,end]);
width_Gabor = zeros(nCh,2);
sd_Gaussian = zeros(nCh,1);
for i=1:nCh
    
    width_Gabor(i,:) = [tbl.fit{i}.Gabor.width80.x1,...
        tbl.fit{i}.Gabor.width80.x2];
    
    sd_Gaussian(i) = tbl.fit{i}.Gaussian.param(2);
 
    [r,c] = ind2sub([ncol,nrow-1],i);
    ha_ch(r,c) = subplot(nrow,ncol,i);
    
    x = tbl.fit{i}.dat.x;
    y = tbl.fit{i}.dat.y;
    se = tbl.fit{i}.dat.se;
    
    errorbar(x,y,se,'LineStyle','none','Marker','o','Color','k');
    
    yy = gabor(xx,tbl.fit{i}.Gabor.param);
    line(xx,yy,'Color',lcolor(1,:),'LineWidth',lineWidth,...
        'DisplayName','Gabor');
    
    str = sprintf('Ch %d',i);
    title(str,'FontSize',14)

    ylim = get(gca,'YLim');
    ylim(1) = 0;
    ylim(2) = ylim(2)*1.1;
    str = sprintf('%3.1f (%3.2f)',tbl.peakX_Gabor(i),tbl.Rsqr_Gabor(i));
    text(xlim(1)+diff(xlim)*0.05,ylim(2)-diff(ylim)*0.05,...
        str,'Color',lcolor(1,:),'FontSize',12);
 
    set(gca,'Box','off','TickDir','out','XLim',xlim,'YLim',ylim)
end

ha_ch = ha_ch';

set(gcf,'CurrentAxes',ha_ch(nrow-1,1))
xlabel('X (deg)','FontSize',14)
ylabel('Firing Rate (spk/s)','FontSize',14)

for i=1:size(ha_ch,1)
    for j=1:size(ha_ch,2)
        ha_ch(i,j).Position(2) = ha_ch(i,j).Position(2) + 0.03;
    end
end

mrg = 1 - (ha_ch(1,1).Position(2)+ha_ch(1,1).Position(4));
yoffset = (ha_ch(1,1).Position(2)-...
    (ha_ch(2,1).Position(2)+ha_ch(2,1).Position(4)))*1.2;
pos = zeros(2,4);
pos(1,1) = ha_ch(1,1).Position(1);
pos(2,1) = ha_ch(1,3).Position(1);
pos(:,3) = ha_ch(1,2).Position(1)+ha_ch(1,2).Position(3) - ...
    ha_ch(1,1).Position(1);
pos(:,2) = mrg;
pos(:,4) = ha_ch(nrow-1,1).Position(2) - yoffset - mrg;

ha = gobjects(1,2);
ha(1) = axes('Position',pos(1,:));
x = 1:nCh;
y = tbl.peakX_Gabor;
errorbar(x,y,y-width_Gabor(:,1),width_Gabor(:,2)-y,'Marker','none',...
    'Color',lcolor(1,:))
idx = tbl.Rsqr_Gabor >= thr;
line(x(idx),y(idx),'Marker','o','MarkerFaceColor',lcolor(1,:),...
    'MarkerEdgeColor',lcolor(1,:),'MarkerSize',8,'LineStyle','none')
line(x(~idx),y(~idx),'Marker','o','MarkerFaceColor','none',...
    'MarkerEdgeColor',lcolor(1,:),'MarkerSize',8,'LineStyle','none')
xlabel('Channel','FontSize',14)

ylabel('X at Peak (deg)','FontSize',14)
title('Fit with Gabor','FontSize',14)

ha(2) = axes('Position',pos(2,:));
y = tbl.peakX_Gaussian;
errorbar(x,y,sd_Gaussian,'Marker','none','Color',lcolor(2,:))
idx = tbl.Rsqr_Gaussian >= thr;
line(x(idx),y(idx),'Marker','o','MarkerFaceColor',lcolor(2,:),...
    'MarkerEdgeColor',lcolor(2,:),'MarkerSize',8,'LineStyle','none')
line(x(~idx),y(~idx),'Marker','o','MarkerFaceColor','none',...
    'MarkerEdgeColor',lcolor(2,:),'MarkerSize',8,'LineStyle','none')
ylim = cell2mat(get(ha,'YLim'));
ylim = [min(ylim(:,1)),max(ylim(:,2))];
xlabel('Channel','FontSize',14)
ylabel(sprintf('%s at Peak',upper(str)),'FontSize',14)
title('Fit with Gaussian','FontSize',14)

set(ha,'Box','off','TickDir','out','XLim',[0.5,nCh+0.5],...
    'YLim',ylim)

end


function G = gabor(x,p)

if isfield(p,'xpower')
    xt = sign(x) .* abs(x).^p.xpower;
    xt0 = sign(p.x0)*abs(p.x0)^p.xpower;
    s = cos(2*pi*p.w*((xt- xt0))+deg2rad(p.phi));
    g = exp(-1*(xt-xt0).^2/(2*p.sigma^2));
else
    s = cos(2*pi*p.w*((x - p.x0))+deg2rad(p.phi));
    g = exp(-1*(x-p.x0).^2/(2*p.sigma^2));
end
G= p.gain * (s .* g) + p.offset;     % values of g range  from -1 to 1,

end



function R_prime = gaussian(param,x)
alpha = param(1);
beta = param(2);
amp = param(3);
baseline = param(4);

r = exp(-1*(x-alpha).^2/beta^2);

R_prime = baseline + amp * (r - min(r));

end