function tbl = fitOneEx(ex)


if ischar(ex)
    load(ex,'ex')
end

% cluster, peakX, width, Rsqr, dat, fit
s = fitEx(ex);
ns = length(s);
validIdx = true(ns,1);
for i=1:length(s)
    if isempty(s{i}) || isempty(s{i}.fit)
        validIdx(i) = false;
    end
end
s = s(validIdx);
ns = length(s);

tbl = [table(string(repmat(s{1}.cluster,ns,1)),'VariableNames',{'cluster'}),...
    array2table(zeros(ns,6),'VariableNames',{'peakX_Gabor','Rsqr_Gabor',...
    'width_Gabor','peakX_Gaussian','Rsqr_Gaussian','sd_Gaussian'}),...
    cell2table(cell(ns,2),'VariableNames',{'dat','fit'})];

validIdx = true(ns,1);
for i=1:length(s)
    tbl.cluster(i) = s{i}.cluster;
    if isempty(s{i}.fit)
        validIdx(i) = false;
    else
        tbl.peakX_Gabor(i) = s{i}.fit.Gabor.max.x;
        tbl.Rsqr_Gabor(i) = s{i}.fit.Gabor.Rsqr;
        tbl.width_Gabor(i) = s{i}.fit.Gabor.width80.width;
        tbl.peakX_Gaussian(i) = s{i}.fit.Gaussian.param(1);
        tbl.Rsqr_Gaussian(i) = s{i}.fit.Gaussian.R_sqr;
        tbl.sd_Gaussian(i) = s{i}.fit.Gaussian.param(2);
        tbl.dat{i} = s{i}.dat;
        tbl.fit{i}.dat = s{i}.fit.dat;
        tbl.fit{i}.Gabor = s{i}.fit.Gabor;
        tbl.fit{i}.Gaussian = s{i}.fit.Gaussian;
    end
end
tbl = tbl(validIdx,:);

end


function output = fitEx(ex)
minN = 3;

stimVar = ex.exp.e1.type;
if ~strcmpi(stimVar,{'x0','y0'})
    error('%s: invalid exp type',stimVar)
end


flds = fieldnames(ex.Trials);
idx = false(length(flds),1);
for i=1:length(flds)
    if isnumeric(ex.Trials(1).(flds{i})) && ...
            length(ex.Trials(1).(flds{i})) == 1
        idx(i) = true;
    end
end

% assign NaN to empty elements in oRate
idx2kill = ~idx;
oRate = struct2table([ex.Trials.oRate]);
cls = oRate.Properties.VariableNames;
oRate = table2cell(oRate);
Trials = [struct2table(rmfield(ex.Trials,flds(idx2kill))),...
    cell2table(oRate,'VariableNames',cls)];


% find stimulus values
idx = Trials.st ~= 0;
x = unique(Trials.(stimVar)(idx));
nx = length(x);

if nx < 5
    output = [];
    return
end

if sum(~idx) > 0
    blankExist = true;
    ns = nx + 1;
else
    blankExist = false;
    ns = nx;
end

ncls = length(cls);
dat = cell(ncls,1);
varNames = {'x','n','meanRate','stRate'};
output = cell(ncls,1);
for j=1:ncls
    mtx = NaN(ns,length(varNames));
    idxNumeric = ~isnan(Trials.(cls{j}));
    for i=1:nx
        idx = idxNumeric & Trials.(stimVar) == x(i);
        mtx(i,:) = [x(i), sum(idx), ...
            mean(Trials.(cls{j})(idx)),...
            std(Trials.(cls{j})(idx))];
    end
    
    if blankExist
        idx = idxNumeric & Trials.st == 0;
        mtx(i+1,:) = [Inf, sum(idx), ...
            mean(Trials.(cls{j})(idx)),...
            std(Trials.(cls{j})(idx))];
    end
    
    dat{j} = array2table([mtx,zeros(ns,1)],'VariableNames',[varNames,'seRate']);
    dat{j}.seRate = dat{j}.stRate ./ sqrt(dat{j}.n);
    
    % fit here
    if std(unique(dat{j}.meanRate)) > 0
        
        if min(dat{j}.n) >= minN
            s1 = fitRF(dat{j},'plot',false);
            s2 = fitGaussian(dat{j},'plot_opt',false);
            output{j}.cluster = cls{j};
            output{j}.dat = dat{j};
            output{j}.fit.dat = s1.dat;
            output{j}.fit.Gabor = s1.fit;
            output{j}.fit.Gaussian = s2;
        end
    end
end

end



function  s = fitRF(tbl,varargin)

idx = ~isinf(tbl.x);
x = tbl.x(idx);
y = tbl.meanRate(idx);
se = tbl.seRate(idx);

meanBlank = [];
if sum(~idx) == 1
    meanBlank = tbl.meanRate(~idx);
    seBlank = tbl.seRate(~idx);
end

plot_opt = false;
opt_fun = 'lsqnonlin';
par_input = struct('w',[],'phi',[],'x0',[],'sigma',[],...
    'gain',[],'offset',[],'xpower',[]);
par_names = fieldnames(par_input);
opts = [];
xtransform = false;
pw = 1;
lim_pwr = 0.8;
res_minmax = 0.0001;
argin = reshape(varargin,2,length(varargin)/2)';
i = contains(argin(:,1),'init','IgnoreCase',true);
if sum(i) == 1
    par_input = argin{i,2};
    argin(i,:) = [];
elseif sum(i) > 1
    error('multiple init')
else
end

for i=1:size(argin,1)
    par_name = lower(argin{i,1});
    par_val = argin{i,2};
    switch par_name
        case 'plot'
            plot_opt = par_val;
        case 'opt_fun'
            opt_fun = par_val;
        case 'w'
            par_input.w = par_val;
        case 'phi'
            par_input.phi = par_val;
        case 'x0'
            par_input.x0 = par_val;
        case 'sigma'
            par_input.sigma = par_val;
        case 'gain'
            par_input.gain = par_val;
        case 'offset'
            par_input.offset = par_val;
        case 'opts'
            opts = par_val;
        case 'xpower'
            pw = par_val;
        case 'xtransform'
            xtransform = par_val;
        case 'lim_pwr'
            lim_pwr = par_val;
        otherwise
    end
end

dat.x = x; dat.y = y; dat.se = se;

% w, phi, x0, sigma, gain ,offset
param_init = [0, 0, 0, 1, 1, 0, pw];
[~, idx] = max(dat.y);
param_init(1) = 1.0e-2;
param_init(3) = dat.x(idx);                         % x0
param_init(4) = (max(dat.x) - min(dat.x)) * 0.5;    % sigma
param_init(5) = max(dat.y) - min(dat.y);            % gain
param_init(6) = mean(dat.y);                        % offset

for i=1:length(par_names)
    if ~isempty(par_input.(par_names{i}))
        param_init(i) = par_input.(par_names{i});
    end
end

y_prime = NaN(size(dat.y));

if strcmpi(opt_fun,'lsqnonlin')
    if isempty(opts)
        opts = optimset('TolX',1e-12,'TolFun',1e-12,'MaxIter',10000,...
            'Display','off','LargeScale','on');
    end
    
    y_range = max(dat.y) - min(dat.y); 
    Au = y_range * 2;
    Al = y_range * 0.001;
    Bu = max(dat.y);
    Bl = -1 * Bu;
    lb = [0, -90, min(dat.x), 0.01, Al, Bl, lim_pwr];
    ub = [3, 90, max(dat.x), max(dat.x) - min(dat.x), Au, Bu, 1];
    
    % limit to 2 cycles/range
    ub(1) = 2/diff(x([1,end]));

    if xtransform
        [param, sse, res, exitflag, output] = ...
            lsqnonlin(@gabor_lsqnonlin2, param_init, lb, ub, opts);
    else
        %param_init(end) = []; lb(end) = []; ub(end) = [];
        [param, sse, res, exitflag, output] = ...
            lsqnonlin(@gabor_lsqnonlin, param_init, lb, ub, opts);
    end
    
elseif strcmpi(opt_fun,'fminsearch')
    if isempty(opts)
        opts = optimset('TolX',1e-12,'TolFun',1e-12,'MaxIter',10000,...
            'Display','off','MaxFunEvals',5000);
    end
    [param, sse, exitflag, output] = ...
        fminsearch(@gabor_fminsearch, param_init(1:6), opts);
    res = y - y_prime;
else
    
end

ssTotal = sum((y - mean(y)).^2);
Rsqr = 1 - sse/ssTotal;

s.dat = dat;

for i=1:length(param)
    s.fit.param.(par_names{i}) = param(i);
end


% find max and min
xx = min(dat.x):res_minmax:max(dat.x);

s.fit.minmax_res = res_minmax;

yy=gabor(xx,s.fit.param);

[r,idx] = max(yy);
s.fit.max.x = xx(idx);
s.fit.max.y = r;


[r,idx] = min(yy);
s.fit.min.x = xx(idx);
s.fit.min.y = r;

options = optimset('TolX',1e-8);
fun = @(x)-1*gabor(x,s.fit.param);
[x_max,fval] = fminbnd(fun,x(1),x(end),options);
if -1*fval > s.fit.max.y
    s.fit.max.x = x_max;
    s.fit.max.y = -1*fval;
end

fun= @(x)gabor(x,s.fit.param);
[x_min,fval] = fminbnd(fun,x(1),x(end),options);
if fval < s.fit.min.y
    s.fit.min.x = x_min;
    s.fit.min.y = fval;
end

s.fit.max.endpoint = false;
s.fit.min.endpoint = false;
if min(abs(x([1,end])-s.fit.max.x)) < res_minmax
    s.fit.max.endpoint = true;
end
if min(abs(x([1,end])-s.fit.min.x)) < res_minmax
    s.fit.min.endpoint = true;
end

% get RF width here
height80 = s.fit.max.y - (s.fit.max.y - s.fit.param.offset) * 0.8;

keepSearching = true;
while keepSearching
    idx_left = find(yy < height80 & xx < s.fit.max.x,1,'last')+1;
    if isempty(idx_left) % left-most end point
        idx = xx < s.fit.max.x;
        if any(diff(yy(idx)) < 0)
            keepSearching = false;
            idx_left = NaN;
        end
        if (s.fit.max.x - xx(1)) > 90
            idx_left = NaN;
            keepSearching = false;
        end
    else
        if idx_left > 1
            keepSearching = false;
        else
            if (s.fit.max.x - xx(idx_left)) > 90
                idx_left = NaN;
                keepSearching = false;
            end
        end
    end
    if keepSearching
        % expand search range by 10 degrees
        xx = xx(1)-10:res_minmax:xx(end);
        yy = gabor(xx,s.fit.param);
    end
end


keepSearching = true;
while keepSearching
    len_xx = length(xx);
    idx_right = find(yy < height80 & xx > s.fit.max.x,1,'first')-1;
    if isempty(idx_right)   % right-most end point
        idx = xx > s.fit.max.x;
        if any(diff(yy(idx))>0)
            keepSearching = false;
            idx_right = NaN;
        end
        if (xx(end) - s.fit.max.x) > 90
            idx_right = NaN;
            keepSearching = false;
        end
    else
        if idx_right < len_xx
            keepSearching = false;
        else
            if (xx(end) - s.fit.max.x) > 90
                idx_right = NaN;
                keepSearching = false;
            end
        end
    end
    if keepSearching
        %xx = min(xx):res_minmax:min([90,max(xx)+(max(xx)-min(xx))]);
        xx = xx(1):res_minmax:xx(end)+10;
        yy = gabor(xx,s.fit.param);
    end
end

if isnan(idx_left)
    s.fit.width80.x1 = NaN;
    s.fit.width80.y1 = NaN;
else
    s.fit.width80.x1 = xx(idx_left);
    s.fit.width80.y1 = yy(idx_left);
end
if isnan(idx_right)
    s.fit.width80.x2 = NaN;
    s.fit.width80.y2 = NaN;
else
    s.fit.width80.x2 = xx(idx_right);
    s.fit.width80.y2 = yy(idx_right);
end
s.fit.width80.width = s.fit.width80.x2 - s.fit.width80.x1;

s.fit.freeParams = par_names;
s.fit.optFun = opt_fun;
s.fit.sse = sse;
s.fit.Rsqr = Rsqr;
s.fit.residual = res;
s.fit.exitflag = exitflag;
s.fit.output = output;
s.fit.g = y_prime;

if plot_opt
    clf
    mrgn = mean(diff(x))*0.2;
    xlim = [xx(1)-mrgn,xx(end)+mrgn];
    line(xx,yy)
    line([s.fit.width80.x1,s.fit.width80.x2],...
        [s.fit.width80.y1,s.fit.width80.y2],'Marker','x',...
        'LineStyle','none','Color','r')
    line(s.fit.min.x,s.fit.min.y,'Marker','*','Color','r')
    line(s.fit.max.x,s.fit.max.y,'Marker','*','Color','r')
    set(gca,'NextPlot','add')
    errorbar(dat.x,dat.y,dat.se,'Marker','o','LineStyle','none','Color','k');
    if ~isempty(meanBlank)
        errorbar(mean(xlim),meanBlank,seBlank,...
            'Marker','o','LineStyle','none','Color',ones(1,3)*0.5);
    end
    set(gca,'TickDir','out','Box','off','XLim',xlim);
end

    function err = gabor_lsqnonlin(param)
        % w, phi, x0, sigma, gain ,offset
        w = param(1);
        phi = param(2);
        x0 = param(3);
        sigma = param(4);
        A = param(5);
        B = param(6);
        
        c = cos(2*pi*w*((x - x0))+deg2rad(phi));
        g = exp(-1*(x-x0).^2/(2*sigma^2));
        y_prime = A * (c .* g) + B;     % values of g range  from -1 to 1,

        err = y - y_prime;
    end

    function err = gabor_lsqnonlin2(param)
        % w, phi, x0, sigma, gain ,offset, x power
        w = param(1);
        phi = param(2);
        x0 = param(3);
        sigma = param(4);
        A = param(5);
        B = param(6);
        pwr = param(7);
        
        xt = sign(x) .* abs(x).^pwr;
        xt0 = sign(x0)*abs(x0)^pwr;
        c = cos(2*pi*w*((xt - xt0))+deg2rad(phi));
        g = exp(-1*(xt-xt0).^2/(2*sigma^2));
        y_prime = A * (c .* g) + B;     % values of g range  from -1 to 1,

        err = y - y_prime;
    end

    function err = gabor_fminsearch(param)
        % w, phi, x0, sigma, gain ,offset
        w = param(1);
        phi = param(2);
        x0 = param(3);
        sigma = param(4);
        A = param(5);
        B = param(6);

        c = cos(2*pi*w*((x - x0))+deg2rad(phi));
        g = exp(-1*(x-x0).^2/(2*sigma^2));
        y_prime = A * (c .* g) + B;     % values of g range  from -1 to 1,
        
        err = sum((y-y_prime).^2);
    end


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



function  s = fitGaussian(tbl,varargin)

plot_opt = true;
param_init = [];
for i=1:length(varargin)/2
    par_name = lower(varargin{(i-1)*2+1});
    par_val = lower(varargin{i*2});
    switch par_name
        case 'plot_opt'
            plot_opt = par_val;
        case 'init'
            param_init = par_val;
        otherwise
    end
end


opts = optimoptions('lsqnonlin','MaxIter',100,...
    'Display','off');

if plot_opt == true
    opts = optimoptions(opts,'OutputFcn',@showfit);
end

idx = ~isinf(tbl.x);
x = tbl.x(idx);
y = tbl.meanRate(idx);

minR = min(y);
[maxR, idx] = max(y);
if isempty(param_init)
    param_init = zeros(1,4);
    param_init(1) = x(idx);
    param_init(2) = diff(x([1,end]))/2;
    param_init(3) = maxR - minR;
    param_init(4) = minR;
end

lb = [-90, 0.1, 0, 0];
ub = [90, 90, param_init(3)*2, maxR];

try
    [param, sse, res, exitflag, output] = ...
        lsqnonlin(@errorfcn, param_init, lb, ub, opts);
catch ME
    s = ME;
    return
end

R_prime = gaussian(param,x);
R_sqr=1-sum((R_prime-y).^2)/sum((y-mean(y)).^2);

s.param = param;
s.R_sqr = R_sqr;
s.fun = 'lsqnonlin';
s.resnorm = sse;
s.residual = res;
s.exitflag = exitflag;
s.output = output;

    function res = errorfcn(param)
        res = y - gaussian(param,x);
    end
    
    function stop = showfit(param, optimValues, state)        
        stop = false;
        
        switch state
            case 'init'
                clf
                ha = axes('Position',[0.2000    0.2000    0.6000    0.6000],...
                    'TickDir','out','Box','off','FontSize',8,...
                    'TickLength',[0.02, 0.025],'Tag','main axis');
                xlabel('X (deg)','FontSize',10)
                ylabel('Firing Rate (spk/s)','FontSize',10)
                
                line(x,y,'Marker','o','MarkerFaceColor','k',...
                    'MarkerEdgeColor','k','LineStyle','none');
                xx = (0:0.01:1)*(max(x)-min(x))+min(x);
                yy = gaussian(param,xx);
                hl=line(xx, yy);
                set(hl,'Tag','fitline');
                ylim = get(ha,'YLim');
                set(ha,'YLim',[0,ylim(2)]);
                drawnow;
            case 'iter'
                ha = findobj('tag','main axis');
                set(ha,'YLimMode','auto')
                hl = findobj('Tag','fitline');
                xx = get(hl,'XData');
                yy = gaussian(param,xx);
                set(hl, 'YData', yy);
                drawnow;
            case 'done'
                ha = findobj('tag','main axis');
                set(ha,'YLimMode','auto')
                stop = true;
            otherwise
        end
    end
end

function R_prime = gaussian(param,x)
alpha = param(1);
beta = param(2);
amp = param(3);
baseline = param(4);

r = exp(-1*(x-alpha).^2/beta^2);

R_prime = baseline + amp * (r - min(r));

end

