clear all;
close all;
clc;

%{
'sample x strains' matrix with relative abundance normalized by total reads
%}
%% LOAD DATA: UNDEFINED FMT
T = readtable('1_UndefinedHuman_20200722_FMT1_MIDAS.xlsx');
samples = string(T.sample_id);
strains = string(T.species_id);
reads = T.count_reads;

[usamples,~,o2u_samples] = unique(samples);
[ustrains,~,o2u_strains] = unique(strains);

N2 = length(ustrains);
M2 = length(usamples);
relabus2 = zeros(N2,M2);
for i = 1:N2
    inds1 = o2u_strains==i;
    for j = 1:M2
        inds2 = o2u_samples==j;
        relabus2(i,j) = sum(reads(inds1&inds2));
    end
end
tot_reads2 = sum(relabus2);
relabus2 = relabus2./sum(relabus2);


%{
Parse sample names to identify the 'mouse' and 'weeks'
only consider samples that:
  - contain "FMT1-W"
  - does NOT contain "PT31"
  - does NOT contain "Cecum"
split at "-", parse week = 2 element, mouse = 3 element

%}
mouses = zeros(M2,1);
weeks = zeros(M2,1);
for i = 1:M2
    sthis = usamples{i};
    if contains(sthis,'FMT1-W')...
            && ~contains(sthis,'Cecum') && ~contains(sthis,'PT31')
        temp = strsplit(sthis,'-');
        weeks(i) = str2num(temp{2}(2:end));
        mouses(i) = str2num(temp{3}(6:end));
    end
end

expts2 = usamples;
mouses2 = mouses;
weeks2 = weeks;
species2 = ustrains;

%% LOAD DATA: BACKFILL1
T = readtable('mbfv1_dataframe_minRelAbund_0.csv');
samples = string(T.sample_id);
strains = string(T.name);
reads = T.relative_abundance;

[usamples,~,o2u_samples] = unique(samples);
[ustrains,~,o2u_strains] = unique(strains);

N = length(ustrains);
M = length(usamples);
relabus = zeros(N,M);
for i = 1:N
    inds1 = o2u_strains==i;
    for j = 1:M
        inds2 = o2u_samples==j;
        relabus(i,j) = sum(reads(inds1&inds2));
    end
end
relabus = relabus./sum(relabus);

weeks = zeros(M,1);
mouses = zeros(M,1);
for i = 1:M
    sample_this = usamples{i};
    if ~contains(sample_this,'Com')...
            && ~contains(sample_this,'GF')...
            && ~contains(sample_this,'OB')...
            && ~contains(sample_this,'Pat')...
            && ~contains(sample_this,'MC')...
            && ~contains(sample_this,'MDC')
        weeks(i) = str2num(sample_this(2));
        mouses(i) = str2num(sample_this(4:strfind(sample_this,'_')-1));
    end
end
species = ustrains;

%{
find and remove all species whose relative abundance is 0 across all samples
%}
inds = find(all(relabus==0,2));
relabus(inds,:) = [];
species(inds) = [];
[N,M] = size(relabus);

expts1 = usamples;
N1 = N;
M1 = M;
weeks1 = weeks;
mouses1 = mouses;
relabus1 = relabus;
species1 = species;

%% COMBINE DATASETS
%-Find union set of species
species1 = species1;
species2 = species2;
[species_all,ia,ib] = union(species1,species2);

N = length(species_all);
M1 = M;
M = M1 + M2;
relabus_all = zeros(N,M);
for i = 1:N
    sthis = species_all{i};
    for j = 1:M1
        ind = find(species1==sthis);
        if ~isempty(ind)
            relabus_all(i,j) = relabus1(ind,j);
        end
    end
    for j = M1+1:M
        ind = find(species2==sthis);
        if ~isempty(ind)
            relabus_all(i,j) = relabus2(ind,j-M1);
        end
    end
end

mousest = [mouses1;mouses2];
weekst = [weeks1;weeks2];
exptst = [expts1;expts2];
Mt = M;

%% LOAD DATA: BACKFILL2
T = readtable('8_backfill2_mbfv2_dataframe_minRelAbund_0_Sunit_backfill2_MIDAS.csv');
samples = string(T.sample_id);
strains = string(T.name);
reads = T.count_reads;

[usamples,~,o2u_samples] = unique(samples);
[ustrains,~,o2u_strains] = unique(strains);

N = length(ustrains);
M = length(usamples);
relabus = zeros(N,M);
for i = 1:N
    inds1 = o2u_strains==i;
    for j = 1:M
        inds2 = o2u_samples==j;
        relabus(i,j) = sum(reads(inds1&inds2));
    end
end
relabus = relabus./sum(relabus);

weeks = zeros(M,1);
mouses = zeros(M,1);
for i = 1:M
    sample_this = usamples{i};
    if contains(sample_this,'SCV2')...
            && ~contains(sample_this,'Cecum')...
            && ~contains(sample_this,'Small')
        temp = strsplit(sample_this,'-');
        weeks(i) = str2num(temp{2}(5:end));
        mouses(i) = str2num(temp{3}(6:end));
    end
end
species = ustrains;

inds = find(all(relabus==0,2));
relabus(inds,:) = [];
species(inds) = [];
[N,M] = size(relabus);

expts3 = usamples;
N3 = N;
M3 = M;
weeks3 = weeks;
mouses3 = mouses;
relabus3 = relabus;
species3 = species;

%% COMBINE DATASETS
%-Find union set of species
species1 = species_all;
species2 = species3;
[species_all,ia,ib] = union(species1,species2);

N = length(species_all);
M = Mt + M3;
relabus1 = relabus_all;
relabus2 = relabus3;
relabus_all = zeros(N,M);
for i = 1:N
    sthis = species_all{i};
    for j = 1:Mt
        ind = find(species1==sthis);
        if ~isempty(ind)
            relabus_all(i,j) = relabus1(ind,j);
        end
    end
    for j = Mt+1:M
        ind = find(species2==sthis);
        if ~isempty(ind)
            relabus_all(i,j) = relabus2(ind,j-Mt);
        end
    end
end

mouses = [mousest;mouses3];
weeks = [weekst;weeks3];
expts = [exptst;expts3];

%% LOAD TAXONOMY
names = readtable('midas_genome_taxonomy.csv');
species_list = names.species_id;
fams_list = names.phylum;
species_list = string(species_list);
fams_list = string(fams_list);
tax_level = 2;

fams = cell(N,1);
for i = 1:N
    sthis = species_all{i};
    inds = find(species_list==sthis,1);
    if ~isempty(inds)
        fams{i} = fams_list{inds};
    end
end
fams = string(fams);
for i = 1:N
    if strcmp(fams{i},'delta/epsilon subdivisions')
        fams{i} = 'Proteobacteria';
    end
end
fams{41} = 'Bacteroidetes';
[ufams,~,o2u_fams] = unique(fams);
Nfams = length(ufams);

ufams_color = zeros(Nfams,3);
for i = 1:Nfams
    sthis = ufams{i};
    if ~isempty(sthis)
        if tax_level==2
            if strcmp(sthis,'Actinobacteria')
                ufams_color(i,:) = [57,83,165]./255;
            elseif strcmp(sthis,'Bacteroidetes')
                ufams_color(i,:) = [13,125,63]./255;
            elseif strcmp(sthis,'Firmicutes')
                ufams_color(i,:) = [238,32,37]./255;
            elseif strcmp(sthis,'Proteobacteria')
                ufams_color(i,:) = [112,48,160]./255;
            elseif strcmp(sthis,'Verrucomicrobia')
                ufams_color(i,:) = [246,134,32]./255;
            else
                ufams_color(i,:) = [200,200,200]./255;
            end
        end
    end
end

%% (Fig 3C) hCom's vs FMTs
%-Pick out data to plot
expt_inds1 = [find(ismember(mouses1,7:9) & weeks1==4)'];
relabus1 = relabus_all(:,expt_inds1);
relabus1 = nanmean(relabus1,2);
expt_inds2 = [find(ismember(mouses2,5:8) & weeks2==4)']+M1;
relabus2 = relabus_all(:,expt_inds2);
relabus2 = nanmean(relabus2,2);
expt_inds3 = [find(ismember(mouses2,9:12) & weeks2==4)']+M1;
relabus3 = relabus_all(:,expt_inds3);
relabus3 = nanmean(relabus3,2);
expt_inds4 = [find(ismember(mouses2,13:16) & weeks2==4)']+M1;
relabus4 = relabus_all(:,expt_inds4);
relabus4 = nanmean(relabus4,2);
expt_inds5 = [find(ismember(mouses3,16:20) & weeks3==4)']+Mt;
relabus5 = relabus_all(:,expt_inds5);
relabus5 = nanmean(relabus5,2);

xthis = 1:5;
relabus_this = [relabus1,relabus5,relabus2,relabus3,relabus4];
Mthis = length(xthis);

%-Set 0 values to minimum value on yaxis
xthis = 1:Mthis;
ymin = 1e-6;
relabus_this(relabus_this<ymin) = ymin;

%-Group by existence before week 4 and add jitter
xwidth = 0.9;
yres = 0.2;

xjitted = zeros(N,Mthis);
for i = 1:Mthis
    ythis = log10(relabus_this(:,i));
    xjits = AddJitterX(ythis,yres,xwidth);
    xjitted(:,i) = xjits + xthis(i);
end

%-Plot
figthis = figure;
hold on;
box on;
ax = gca;
ax.ActivePositionProperty = 'position';

%-Plot threshold
xtemp = [-10,10];
ytemp = [1,1].*1e-4;
plot(xtemp,ytemp,':','color',[0.2,0.2,0.2,0.5],'linewidth',0.5);

%-Plot species
for i = 1:N
    xtemp = xjitted(i,:);
    ytemp = relabus_this(i,:);
    inds = find(ytemp==ymin);
    xtemp(inds) = [];
    ytemp(inds) = [];
    scatter(xtemp,ytemp,30,...
        ufams_color(o2u_fams(i),:),'filled',...
        'markerfacealpha',1,'markeredgecolor','k',...
        'linewidth',0.25);
end

set(gca,'yscale','log');
ylim([ymin,1]);
yticks([1e-6,1e-4,1e-2,1e0]);
ylabel('Relative abundance');
xticks(xthis);
xticklabels({'hCom1','hCom2','Hum1','Hum2','Hum3'});
xlim([0.5,5.5]);

set(gca,'TickDir','out');
set(gca,'linewidth',0.75,'ticklength',[0.0145,0.00]);
set(gca,'layer','top');
set(gcf,'units','points','InnerPosition',[100,100,220,180]);
set(gca,'fontsize',10);
xl = get(gca,'XLabel');
set(xl,'FontSize',10);
yl = get(gca,'YLabel');
set(yl,'FontSize',10);
set(gca,'FontName','Arial')
set(gca,'xminortick','off','yminortick','off');
ax.XAxis.TickLength = [0,0];

set(figthis,'Renderer','painter');
% saveas(figthis,'FigComps\3d.svg','svg');