% this script load data of (x,y) position of colonies into MATLAB and
% calculate their associated summary statistics 

clear
clc
close all

load("samples v1 i=50 23-Jun-2025 14:08:45.mat")
samples1 = samples(1:50,:);

samples_combined = samples1;

parameters = vertcat(samples_combined{:,3});
%%
disp("Data Loaded!")
N = length(samples_combined);
Tslices = 1;
outer_area = nan(N,Tslices);
total_area = nan(N,Tslices);
compactness = nan(N,Tslices);
rmin = nan(N,Tslices);
rmax = nan(N,Tslices);
rmean = nan(N,Tslices);
branch_count = nan(N,Tslices);
branch_length = nan(N,Tslices);
store_file_name = cell(N,1);

galaxy_area = 72897;

tic
kk = 1;
for ii = 1%1:N
    xF = samples_combined{ii,1};
    yF = samples_combined{ii,2};
    
    min_cell = 0;
    max_cell = 60000;
    
    Ncells = get_n_cells(min_cell,max_cell,xF,yF,galaxy_area);

    cc = 1;
    
    tt = Ncells;
    x = xF(1:tt,:);
    y = yF(1:tt,:);
    
    Telong = round(samples_combined{ii,3}(1),3);
    p2sProb = round(samples_combined{ii,3}(2),3);
    s2pProb = round(samples_combined{ii,3}(3),3);
    pc = round(samples_combined{ii,3}(4),3);
    pa = round(samples_combined{ii,3}(5),3);
    
    % saving colony as tiff
    file_name = "S"+ii+"_N="+tt+"_Telong="+Telong+"_p2sProb="+p2sProb+"_s2pProb="...
                +s2pProb+"_pc="+pc+"_pa="+pa; 

    % save colony as matrix
    [I_hist,~,~] = histcounts2(x+225, y+300,0:2:550,0:2:600);
    I = imresize(I_hist,2)>0.5;
    I = ~bwareaopen(~I, 10); % fill gaps in colony that has 10 or less pixels
    % 
    % [outer_area(ii,cc), total_area(ii,cc)] = get_colony_area_rmean(I);
    % compactness(ii,cc) = get_compactness(I);
    % [~, ~, rmin(ii,cc), rmax(ii,cc), rmean(ii,cc)] = get_radii_rmean(I);
    % [branch_count(ii,cc), branch_length(ii,cc)] = get_branch_length_rmean(I);
    
    cc = cc+1;
    
    store_file_name{ii} = I<0.5;
    kk = kk+1;
    if mod(ii,2) == 0
        fprintf("%0.2f completed and took %0.2f\n",ii/N*100,toc)
    end

end
toc

%%
figure
idx = [1 3 4 10 11 36 41 42 43];
% idx = [1 3 4 10 11 29 36 41 42 43];
disp_colony = store_file_name(idx);

Iexp = imread("colonies_2_t6.tiff")<0.5;

Iexp2 = zeros(size(Iexp)+5);

Iexp2(3:end-3, 3:end-3) = Iexp;

disp_colony{5} = Iexp2;

montage(disp_colony([2,5,8]),'ThumbnailSize',[500 inf],'Size',[1,3])

% exportgraphics(gcf,'panel_galaxy_06Aug2025_w_exp.png','Resolution',200)


%%
