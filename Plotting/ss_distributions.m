% this script takes the summary statistic of the final time for galaxy
% yeast and show boxplots and normalise the summary statisitics near a unit
% interval.
%
% Kai Li
% 15 July 2025

clear
clc
close all

load("galaxy_yeast_single_09May2025.mat");

load("ss_galaxy_1_12-May-2025 07:32:46.mat")

branch_length0 = branch_length;
compactness0 = compactness;
outer_area0 = outer_area;
rmax0 = rmax;
rmean0 = rmean;
total_area0 = total_area;
parameters0  = parameters;

load("ss_galaxy_1_21-Jun-2025 20:08:15.mat")

branch_length1 = branch_length;
compactness1 = compactness;
outer_area1 = outer_area;
rmax1 = rmax;
rmean1 = rmean;
total_area1 = total_area;
parameters1  = parameters;

load("ss_galaxy_1_23-Jun-2025 06:18:26.mat")

branch_length2 = branch_length;
compactness2 = compactness;
outer_area2 = outer_area;
rmax2 = rmax;
rmean2 = rmean;
total_area2 = total_area;
parameters2  = parameters;

load("ss_galaxy_1_23-Jun-2025 14:25:44.mat")

branch_length3 = branch_length;
compactness3 = compactness;
outer_area3 = outer_area;
rmax3 = rmax;
rmean3 = rmean;
total_area3 = total_area;
parameters3  = parameters;

load("ss_galaxy_1_24-Jun-2025 10:44:47.mat")

branch_length4 = branch_length;
compactness4 = compactness;
outer_area4 = outer_area;
rmax4 = rmax;
rmean4 = rmean;
total_area4 = total_area;
parameters4  = parameters;

clear branch_length compactness outer_area rmax rmean total_area parameters

parameters = [parameters0; parameters1; parameters2; parameters3; parameters4];
branch_length = [branch_length0; branch_length1; branch_length2; branch_length3; branch_length4];
compactness = [compactness0; compactness1; compactness2; compactness3; compactness4];
outer_area = [outer_area0; outer_area1; outer_area2; outer_area3; outer_area4];
rmax = [rmax0; rmax1; rmax2; rmax3; rmax4];
rmean = [rmean0; rmean1; rmean2; rmean3; rmean4];
total_area = [total_area0; total_area1; total_area2; total_area3; total_area4];


%%

df = [branch_length(:,end), compactness(:,end), outer_area(:,end), rmax(:,end),...
      rmean(:,end), total_area(:,end)];

df_exp = [branch_length_galaxy, compactness_galaxy, outer_area_galaxy, rmax_galaxy,...
          rmean_galaxy, total_area_galaxy];

max_val = max(df);

df_norm = df./max_val;
df_exp_norm = df_exp./max_val;

ss_names_combined = ["Branch Count", "Compactness", "Outer Area", "Rmax",...
                "Rmean", "Total Area"];

boxplot(df_norm,'Labels',ss_names_combined)
hold on
plot(df_exp_norm)
% exportgraphics(gcf,'ss_galaxy_single_09May2025_1417.pdf','ContentType','vector')

%%

bw = [0.05 0.05 0.05 0.05 0.05 0.0005];

figure('Position',[10,10,1200,800])
xLab = ["(a)","(b)","(c)","(d)","(e)","(e)","(f)"];

ss_names_hist = ["Branch Count, $I_B$", "Compactness, $I_C$", "Outer Area, $I_{A_F}$",...
                 "Maximum Radius, $I_{R_{Max}}$", "Mean Radius, $I_{R_{Mean}}$",...
                 "Total Area, $I_A$"];

r_idx = {1:150, 1201:1350, 1401:1650, 1651:1900, 1901:2050};

for ii = 1:6
    nexttile
    p = [];
    cc = 1;
    for r = [1 3 5]
        hold on
        % histogram(df_norm(r_idx{r},ii),'EdgeColor','none','BinWidth',bw(ii),'FaceAlpha',0.5)
        [f,xi] = ksdensity(df_norm(r_idx{r},ii),'Support',[-0.01,1]); 

        p{cc} = plot(xi,f);
        
        xline(df_exp_norm(ii),'LineWidth',2);
        set(gca,'FontSize',24)
        xlabel(xLab(ii),'FontSize',36)
        ylabel(ss_names_hist(ii),'FontSize',28)
        if ii == 6 
            xlim([0.992 1])
        end
      
        cc = cc+1;
        box on
    end

    if ii == 1
        legend([p{1},p{2},p{3}],"Initial Round", "Round 2", "Round 4")
    end

end

% exportgraphics(gcf,"galaxy_ss_histogram_24June2025.pdf",'ContentType','vector')


%%

bw = [0.05 0.05 0.05 0.05 0.05 0.0005];


xLab = ["(a)","(b)","(c)","(d)","(e)","(f)"];

ss_names_hist = ["Branch Count, $I_B$", "Compactness, $I_C$", "Outer Area, $I_{A_F}$",...
                 "Maximum Radius, $I_{R_{Max}}$", "Mean Radius, $I_{R_{Mean}}$",...
                 "Total Area, $I_A$"];

r_idx = {1:150, 1201:1350, 1401:1650, 1651:1900, 1901:2050};

col = lines;

cc = 1;

for r = 1:5
    pause(0.1)
    figure('Position',[10,10,1200,800])
    for ii = 1:5
        nexttile
        hold on
        histogram(df_norm(r_idx{r},ii),'EdgeColor','none','FaceAlpha',0.5,'FaceColor',col(cc,:))

        xline(df_exp_norm(ii),'LineWidth',2);
        set(gca,'FontSize',24)
        title(xLab(ii),'FontSize',36)
        xlabel(ss_names_hist(ii),'FontSize',28)
        ylabel("Frequency","FontSize",28)
        if ii == 6 
            xlim([0.992 1])
        else
            xlim([0,1])
        end
        box on
    end
    cc = cc+1;
    % exportgraphics(gcf,"galaxy_ss_histogram_round_"+(cc-2)+".pdf",'ContentType','vector')
end


% 

%% makes csv for simulated data

parameter_name = ["nStar"; "p2sProb"; "s2pProb"; "gamma"; "pa"; "angle_prolif"];
names_combined = [parameter_name; ss_names_combined'];

logit_parameters(:,1:5) = log(parameters(:,1:5)./(1-parameters(:,1:5)));

logit_parameters(:,6) = log( (parameters(:,6)/(pi/20))./(1-(parameters(:,6)/(pi/20))) );

df_final = [logit_parameters, df_norm];

df_table = array2table(df_final,"VariableNames",names_combined);

% writetable(df_table,"ss_galaxy_r3_23June2025.csv")

%% make csv for experimental data

exp_table = array2table(df_exp_norm,'VariableNames',ss_names_combined);
% writetable(exp_table,"ss_galaxy_exp_r3_ss_23June2025.csv")
