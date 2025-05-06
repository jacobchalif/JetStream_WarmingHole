addpath(genpath("Functions/"))

option = 2;
    

mov = 7;
lw = 1.25;
YL = [-3.5 3.5];
titFS = 18;
legFS = 14;
labFS = 16;
baseFS = 12;
        

SOM_ID = '20240318_NCEP_19802010';
load(['SOMs/' SOM_ID '.mat']);
nclasses = DIM_1*DIM_2;

dtasets = {'NCEP','20CRv3','ERA5','ERA20c','JRA55','MERRA2','NCEP2'};
legDta = dtasets;
starts = [1948 1900 1940 1900 1958 1980 1979 1915];
ends = [2023 2015 2023 2010 2023 2023 2023 2018];
cols = colororder();
cols(end+1,:) = [0 0 0];
overlap = [max(starts) min(ends)];
ssn = 1; % 1 = DJF



%% Selecting wavy classes
LWA = readtable("Data/Processed Reanalysis/LWA/NCEP.csv");
MCI = readtable("Data/Processed Reanalysis/MCI300/NCEP.csv");
classes = readmatrix(['Data/Processed Reanalysis/SOM/' SOM_ID '/master.csv']);
N = findDayIndex(2023,12,31,1948);
dates = nan(N,3);
for year = 1948:2023
    l = isLeap(year);
    m = monthDays(l);
    i1 = findDayIndex(year,1,1,1948);
    i2 = findDayIndex(year,12,31,1948);
    years(i1:i2) = year;
    dates(i1:i2,1) = year;
    for i = 1:12
        if i == 1
            dates(i1:(i1+m(i)-1),2) = i;
            dates(i1:(i1+m(i)-1),3) = 1:m(i);
        else
            dates((i1+sum(m(1:(i-1)))):(i1+sum(m(1:(i-1)))+m(i)-1),2) = i;
            dates((i1+sum(m(1:(i-1)))):(i1+sum(m(1:(i-1)))+m(i)-1),3) = 1:m(i);
        end
    end
end
winterI = ((dates(:,1) >= 1980) & (dates(:,1) <= 2010)) & ((dates(:,2) == 12) | (dates(:,2) <= 2));
dates = dates(winterI,:);
LWA = LWA(winterI,:);
MCI = MCI(winterI,:);


classNames = strings(nclasses,1);
c1 = 0;
c2 = 0;
for class = 1:nclasses
    classNames(class) = strcat("(", num2str(c1), ", " , num2str(c2),")");
    if c2 < (DIM_2-1)
        c2 = c2 + 1;
    else
        c2 = 0;
        c1 = c1 + 1;
    end
end

%% Wavy days
threshold = 0.8;
labFS = 18;
baseFS = 14;
figure('Position',[121 343 1100 700]);
ha = tight_subplot(DIM_1,DIM_2,[.05 .01],[.08 .06],[.06 .02]);  % (Nh, Nw, gap, marg_h, marg_w)
for class = 1:nclasses
    classI = classes(:,4) == class;
    n = sum(classI);
    LWA_class = LWA.LWA_Anticylonic(classI) - LWA.LWA_Cylonic(classI);
    MCI_class = MCI.MCI(classI);
    if ((sum(LWA_class<0)/n) >= threshold) && ((sum(MCI_class<0)/n) >= threshold)
        col = 'r';
    else
        col = 'k';
    end
    axes(ha(class))
    scatter(MCI_class,LWA_class/1e6,15,col,'filled')
    xlim([-1 1])
    ylim([-50 50])
    xline(0)
    yline(0)
    set(gca,'FontSize',baseFS,'FontName','Avenir')
    if class >= 29
        xlabel("MCI",'FontSize',labFS)
    else
        xticklabels([])
    end
    if round((class-1)/7) == (class-1)/7
        ylabel("LWA",'FontSize',labFS)
    else
        yticklabels([])
    end
    title(classNames(class),'FontSize',labFS)
end
exportgraphics(gcf,'Figures/Supplement/SI4.png','Resolution',450)

