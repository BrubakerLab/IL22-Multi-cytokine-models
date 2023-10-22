clear
clc
close all force
mkdir Results
%----------------------------------------------------------------------------------------------------%
% Results Summary
%----------------------------------------------------------------------------------------------------%
%{
______________________________________________________________________________________________________
SINGLE-CYTOKINE MODELING (GENERALIZED LINEAR MODEL GLM)
______________________________________________________________________________________________________

mIL22 GLM: Significant regulation by mIL22 of 
- IL8, VEGF, IP-10, IL1-RA
- Borderline (p < 0.10) RANTES, PDGF_BB, IL17, GM_CSF, Eotaxin, IL_1beta

geIL22 GLM: Significant regulation by geIL22 of:
- IL8, VEGF, IP-10, IL1-RA, IL6, IL5, IL17
- Borderline (p < 0.10) IL15 and IL9

Pure/Isolated core IL22 response (IL8, VEGF, IP-10, IL1-RA) conserved between
mIL22 and geIL22 conditions. Expanded list includes IL17, a benefitial
cytokine associated with Th17 cell recruitment to inflamed epithelium.

Though core regulatory response is conserved, specific conditions (IL22,
LPS, LPSxIL22) that significantly associate with each cytokine are not
similar indicating some difference in the regulatory behavior, likely due to off-target effects.  

Off-target effects of geIL22 likely due to secreted microbial products
include IL15, IL9, IL5, and IL6. Some effects may be beneficial, others may
be problematic (need to do literature dive). 

Secondary effects of mIL22 include Eotaxin, IL1b, GMCSF, RANTES, PDGF-BB.
Unclear if these are beneficial or problematic (need to do literature dive)

______________________________________________________________________________________________________
PARTIAL LEAST SQUARES DISCRIMINANT ANALYSIS (MULTI-CYTOKINE MODELS)
______________________________________________________________________________________________________

In both the mIL22 and geIL22 multivariate models, two axes emerge, one
driven by the IL22 stimulus and the other by the LPS stimulus. critically,
conditions group together and the plots across both experiments mirror each
other, indicating similar overall cytokine networks. 

mIL22 PLS-DA significant cytokines (accounting for correlation structure of
all cytokines at once)
- IL8, VEGF, MCP1, IP-10, PDGF-BB, G-CSF, IL1RA
- Need to examine LV3 properties...
- Need to examine LV1 and LV2 individually significant cytokines...

geIL22 PLSDA Significant Cytokiens
- IL17, IP_10, VEGF, MIP_1alfa, MCP_1, IL12, IL_1ra
- Need to examine LV3 properties...
- Need to examine LV1 and LV2 individually significant cytokines...

SUMMARY: Core signals of IP-10, VEGF, MCP-1, IL1RA are conserved in
multivariate modeling of responses to mIL22 and geIL22. 

The potentially problematic off-target effects of geIL22 (IL15, IL9, IL5, and IL6) are
significant on LV1 alone in the geIL22 model, all of which are loaded
toward the geIL-22 conditions rather than the LPS conditions in the scores
plot, indicating these signals are primarily associated with the
therapeutic molecule rather than the inflammatory stimulus. 

In mIL22 models, IL15 was excluded due to low abundance, but IL9, IL5, and
IL6 are also primarily associated with the therapeutic molecule, mIL22, but
not as strongly/significatnly regualted as compared with geIL22. 

The primary advantage that the geIL22 system seems to have over the mIL22
stimulus is the IL-17 induction and its association with the therapeutic
molecule. IL-17 is associated with recruitment of T-helper 17 (Th17) cells to sites of
inflammation in the gut and has a protective effect. 

An interesting experiment/prediction of the cell culture studies is that if
we were to put the microbe and pure IL22 into mouse gut, we would observe
Th17 cell recruitment to site of inflammation in an IBD mouse model and 
production of IL-17 by the mouse epithelium. 

%}

%----------------------------------------------------------------------------------------------------%
% Hyperparameters 
%----------------------------------------------------------------------------------------------------%
resultsHeaders = { 'Model p-value';...
                   'F statistic';
                   'Intercept beta';
                   'IL22 beta';
                   'LPS beta';
                   'IL22 x LPS beta';
                   'Intercept se';
                   'IL22 se';
                   'LPS se';
                   'IL22 x LPS se';
                   'Intercept Tstat';
                   'IL22 Tstat';
                   'LPS Tstat';
                   'IL22 x LPS Tstat';
                   'Intercept p-value';
                   'IL22 p-value';
                   'LPS p-value';
                   'IL22 x LPS p-value';};



%----------------------------------------------------------------------------------------------------%
% Data Import
%----------------------------------------------------------------------------------------------------%

data_mIL22  = readtable('Final_Data/V2_cleanedLuminex_mIL22_experiment.txt','Delimiter','\t','ReadRowNames',1);
data_geIL22 = readtable('Final_Data/V2_cleanedLuminex_geIL22_experiment.txt','Delimiter','\t','ReadRowNames',1);
data_all    = readtable('Final_Data/V3_cleanedLuminex_allData_experiment.txt','Delimiter','\t','ReadRowNames',1);


%----------------------------------------------------------------------------------------------------%
% Processing
%----------------------------------------------------------------------------------------------------%

mIL22_avgs  = mean(data_mIL22{1:3,5:end}+1);
geIL22_avgs = mean(data_geIL22{1:3,5:end}+1);
aIL22_avgs  = mean(data_all{1:3,6:end}+1);

data_mIL22{:,5:end} = log2((data_mIL22{:,5:end}+1)./mIL22_avgs);
data_geIL22{:,5:end} = log2((data_geIL22{:,5:end}+1)./geIL22_avgs);
data_all{:,6:end} = log2((data_all{:,6:end}+1)./aIL22_avgs);


[~,c1] = size(data_mIL22);
[~,c2] = size(data_geIL22);

%----------------------------------------------------------------------------------------------------%
% Conditions_GLM
%----------------------------------------------------------------------------------------------------%
counter = 1;
for i = 5:c1
    mdltbl_mIL22 = data_mIL22(:,[2 4 i]);
    mdl_mIL22 = fitglm(mdltbl_mIL22,'interactions');

    % Create a table of results of mIL22
    mIL22_results(counter,1)       = mdl_mIL22.coefTest; % model p-value
    mIL22_results(counter,2)       = mdl_mIL22.devianceTest{2,3}; % model p-value
    mIL22_results(counter,[3:6])   = mdl_mIL22.Coefficients{:,1}'; % model coefficients
    mIL22_results(counter,[7:10])  = mdl_mIL22.Coefficients{:,2}'; % model coefficient standard deviation
    mIL22_results(counter,[11:14]) = mdl_mIL22.Coefficients{:,3}'; % model coefficient T statistic
    mIL22_results(counter,[15:18]) = mdl_mIL22.Coefficients{:,4}'; % model coefficient p-values
    
    counter = counter +1;

end
mIL22_T = splitvars(table(mIL22_results,'RowNames',data_mIL22.Properties.VariableNames(5:end)));
mIL22_T.Properties.VariableNames = resultsHeaders;
writetable(mIL22_T,'V3_Results/mIL22_glm_ResultsTable.txt','WriteRowNames',1);

counter = 1;
for i = 5:c2
    mdltbl_geIL22 = data_geIL22(:,[2 4 i]);
    mdl_geIL22 = fitglm(mdltbl_geIL22,'interactions');

    % Create a table of results of geIL22
    geIL22_results(counter,1)       = mdl_geIL22.coefTest; % model p-value
    geIL22_results(counter,2)       = mdl_geIL22.devianceTest{2,3}; % model p-value
    geIL22_results(counter,[3:6])   = mdl_geIL22.Coefficients{:,1}'; % model coefficients
    geIL22_results(counter,[7:10])  = mdl_geIL22.Coefficients{:,2}'; % model coefficient standard deviation
    geIL22_results(counter,[11:14]) = mdl_geIL22.Coefficients{:,3}'; % model coefficient T statistic
    geIL22_results(counter,[15:18]) = mdl_geIL22.Coefficients{:,4}'; % model coefficient p-values
    
    counter = counter +1;

end

geIL22_T = splitvars(table(geIL22_results,'RowNames',data_geIL22.Properties.VariableNames(5:end)));
geIL22_T.Properties.VariableNames = resultsHeaders;
writetable(geIL22_T,'V3_Results/geIL22_glm_ResultsTable.txt','WriteRowNames',1);

%----------------------------------------------------------------------------------------------------%
% PLSDA Model: X: Cytokines, Y:[mIL22 LPS]
%----------------------------------------------------------------------------------------------------%

% PLSDA Model, Multi-Y, mIL22 only
Y_mIL22                    = data_mIL22{:,[2 4]};
X_mIL22                    = zscore(data_mIL22{:,5:end});
[~,~,~,~,~,mPCTVAR,mMSE,~] = plsregress(X_mIL22,Y_mIL22,8,'cv',8);;

% Select Number of LV's
figure
plot(1:8,cumsum(100*mPCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');
title('mIL22 Var Explained')
saveas(gcf,'V3_Results/mIL22_choose_nComp_pctVarY','epsc')

figure
plot(1:8,mMSE(2,1:end-1),'-bo');
xlabel('Number of PLS components');
ylabel('MSE Y');
title('mIL22 MSE')
saveas(gcf,'V3_Results/mIL22_choose_nComp_MSE','epsc')

% Train a 3 LV model
[mXL,mYL,mXS,mYS,mBETA,mPCTVAR,mMSE,mSTATS] = plsregress(X_mIL22,Y_mIL22,3);

figure
gscatter(mXS(:,1),mXS(:,2),data_mIL22.Condition,[],'.',50)
title('X Scores, PLSDA Model mIL22')
xlabel([{mPCTVAR(1,1) '%VarExp LV1'}])
ylabel([{mPCTVAR(1,2) '%VarExp LV2'}])
saveas(gcf,'V3_Results/mIL22_plsda_Xscores','epsc')

figure
gscatter(mYS(:,1),mYS(:,2),data_mIL22.Condition,[],'.',50)
title('Y Scores, PLSDA Model mIL22')
xlabel([{mPCTVAR(2,1) '%VarExp LV1'}])
ylabel([{mPCTVAR(2,2) '%VarExp LV2'}])
saveas(gcf,'V3_Results/mIL22_plsda_Yscores','epsc')


mIL22_vipG       = pls_vip(mSTATS,mXS,mXL,mYL,[1:3]);
for i = 1:3
    mIL22_vipsByLV(:,i) = pls_vip(mSTATS, mXS(:,i), mXL(:,i), mYL,i);
end

tbl_mIL22_vip   = splitvars(table(mIL22_vipG,mIL22_vipsByLV,'RowNames',data_mIL22.Properties.VariableNames(5:end)));
tbl_mIL22_vip.Properties.VariableNames = {'Global';'LV1';'LV2';'LV3'};
writetable(tbl_mIL22_vip,'V3_Results/mIL22_vipScores.txt','WriteRowNames',1);

mVipLegend = tbl_mIL22_vip.Global;
mVipLegend(mVipLegend<1) = 0;
mVipLegend(mVipLegend>1) = 1;

figure
gscatter(mXL(:,1),mXL(:,2),mVipLegend,[],'.',50)
title('Cytokine Loadings, PLSDA Model mIL22')
xlabel([{mPCTVAR(2,1) '%VarExp LV1'}])
ylabel([{mPCTVAR(2,2) '%VarExp LV2'}])
xline(0)
yline(0)
text( mXL(:,1),mXL(:,2) , cellstr(tbl_mIL22_vip.Properties.RowNames(:)) )
saveas(gcf,'V3_Results/mIL22_plsda_Loadings','epsc')

[~,iz] = sort(tbl_mIL22_vip{:,1});
figure
barh(tbl_mIL22_vip{iz,1})
yticks([1:height(tbl_mIL22_vip)])
yticklabels(tbl_mIL22_vip.Properties.RowNames(iz))
xline(1)
saveas(gcf,'V3_Results/mIL22_plsda_vipScores','epsc')

%----------------------------------------------------------------------------------------------------%
% PLSDA Model: X: Cytokines, Y:[geIL22 LPS]
%----------------------------------------------------------------------------------------------------%

% PLSDA Model, Multi-Y, geIL22 only
Y_geIL22 = data_geIL22{:,[2 4]};
X_geIL22 = zscore(data_geIL22{:,5:end});

[~,~,~,~,~,gePCTVAR,geMSE,~] = plsregress(X_geIL22,Y_geIL22,8,'cv',8);;

% Select Number of LV's
figure
plot(1:8,cumsum(100*gePCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');
title('geIL22 Var Explained')
saveas(gcf,'V3_Results/geIL22_choose_nComp_pctVarY','epsc')

figure
plot(1:8,geMSE(2,1:end-1),'-bo');
xlabel('Number of PLS components');
ylabel('MSE Y');
title('geIL22 MSE')
saveas(gcf,'V3_Results/geIL22_choose_nComp_MSE','epsc')


% Train a 3 LV model
[geXL,geYL,geXS,geYS,geBETA,gePCTVAR,geMSE,geSTATS] = plsregress(X_geIL22,Y_geIL22,3);

figure
gscatter(geXS(:,1),geXS(:,2),data_geIL22.Condition,[],'.',50)
title('geIL22 X Scores, PLSDA Model')
xlabel([{gePCTVAR(1,1) '%VarExp LV1'}])
ylabel([{gePCTVAR(1,2) '%VarExp LV2'}])
saveas(gcf,'V3_Results/geIL22_plsda_Xscores','epsc')

figure
gscatter(geYS(:,1),geYS(:,2),data_geIL22.Condition,[],'.',50)
title('geIL22 Y Scores, PLSDA Model')
xlabel([{gePCTVAR(2,1) '%VarExp LV1'}])
ylabel([{gePCTVAR(2,2) '%VarExp LV2'}])
saveas(gcf,'V3_Results/geIL22_plsda_Yscores','epsc')


geIL22_vipG       = pls_vip(geSTATS,geXS,geXL,geYL,[1:3]);

for i = 1:3
    geIL22_vipsByLV(:,i) = pls_vip(geSTATS, geXS(:,i), geXL(:,i), geYL,i);
end

tbl_geIL22_vip   = splitvars(table(geIL22_vipG,geIL22_vipsByLV,'RowNames',data_geIL22.Properties.VariableNames(5:end)));
tbl_geIL22_vip.Properties.VariableNames = {'Global';'LV1';'LV2';'LV3'};
%disp(tbl_geIL22_vip)
writetable(tbl_geIL22_vip,'V3_Results/geIL22_vipScores.txt','WriteRowNames',1);

geVipLegend = tbl_geIL22_vip.Global;
geVipLegend(geVipLegend<1) = 0;
geVipLegend(geVipLegend>1) = 1;

figure
gscatter(geXL(:,1),geXL(:,2),geVipLegend,[],'.',50)
title('Cytokine Loadings, PLSDA Model geIL22')
xlabel([{gePCTVAR(2,1) '%VarExp LV1'}])
ylabel([{gePCTVAR(2,2) '%VarExp LV2'}])
xline(0)
yline(0)
text( geXL(:,1),geXL(:,2) , cellstr(tbl_geIL22_vip.Properties.RowNames(:)) )
saveas(gcf,'V3_Results/geIL22_plsda_Loadings','epsc')

[~,iz] = sort(tbl_geIL22_vip{:,1});
figure
barh(tbl_geIL22_vip{iz,1})
yticks([1:height(tbl_geIL22_vip)])
yticklabels(tbl_geIL22_vip.Properties.RowNames(iz))
xline(1)
saveas(gcf,'V3_Results/geIL22_plsda_vipScores','epsc')


%----------------------------------------------------------------------------------------------------%
% PLSDA Model: X: Cytokines, Y:[mIL22 geIL22 LPS]
%----------------------------------------------------------------------------------------------------%

% PLSDA Model, Multi-Y, all conditions only
Y_all                       = data_all{:,[2 3 5]};
X_all                       = zscore(data_all{:,6:end});
[~,~,~,~,~,aPCTVAR,aMSE,~] = plsregress(X_all,Y_all,8,'cv',8);;

% Select Number of LV's
figure
plot(1:8,cumsum(100*aPCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');
title('All Data Var Explained')
saveas(gcf,'V3_Results/allData_choose_nComp_pctVarY','epsc')

figure
plot(1:8,aMSE(2,1:end-1),'-bo');
xlabel('Number of PLS components');
ylabel('MSE Y');
title('All Data MSE')
saveas(gcf,'V3_Results/allData_choose_nComp_MSE','epsc')

% Train a 3 LV model
[aXL,aYL,aXS,aYS,aBETA,aPCTVAR,aMSE,aSTATS] = plsregress(X_all,Y_all,3);

figure
gscatter(aXS(:,1),aXS(:,2),data_all.Condition,[],'.',50)
title('X Scores, PLSDA Model All Data')
xlabel([{aPCTVAR(1,1) '%VarExp LV1'}])
ylabel([{aPCTVAR(1,2) '%VarExp LV2'}])
saveas(gcf,'V3_Results/allData_plsda_Xscores','epsc')

figure
gscatter(aYS(:,1),aYS(:,2),data_all.Condition,[],'.',50)
title('Y Scores, PLSDA Model All Data')
xlabel([{aPCTVAR(2,1) '%VarExp LV1'}])
ylabel([{aPCTVAR(2,2) '%VarExp LV2'}])
saveas(gcf,'V3_Results/allData_plsda_Yscores','epsc')


all_vipG       = pls_vip(aSTATS,aXS,aXL,aYL,[1:3]);
for i = 1:3
    all_vipsByLV(:,i) = pls_vip(aSTATS, aXS(:,i), aXL(:,i), aYL,i);
end

tbl_all_vip   = splitvars(table(all_vipG,all_vipsByLV,'RowNames',data_all.Properties.VariableNames(6:end)));
tbl_all_vip.Properties.VariableNames = {'Global';'LV1';'LV2';'LV3'};
writetable(tbl_all_vip,'V3_Results/allData_vipScores.txt','WriteRowNames',1);

aVipLegend = tbl_all_vip.Global;
aVipLegend(aVipLegend<1) = 0;
aVipLegend(aVipLegend>1) = 1;

figure
gscatter(aXL(:,1),aXL(:,2),aVipLegend,[],'.',50)
title('Cytokine Loadings, PLSDA Model mIL22')
xlabel([{aPCTVAR(2,1) '%VarExp LV1'}])
ylabel([{aPCTVAR(2,2) '%VarExp LV2'}])
xline(0)
yline(0)
text( aXL(:,1),aXL(:,2) , cellstr(tbl_all_vip.Properties.RowNames(:)) )
saveas(gcf,'V3_Results/allData_plsda_Loadings','epsc')






%----------------------------------------------------------------------------------------------------%
% Specific bar plots of cytokines
%{
mIL22 PLS-DA significant cytokines:
- IL8, VEGF, MCP1, IP-10, PDGF-BB, G-CSF, IL1RA

%}
%----------------------------------------------------------------------------------------------------%


% mIL22 IL8
idx = 7;
disp(data_mIL22.Properties.VariableNames(idx))
scale = mean(data_mIL22{1:3,idx});
plotData = {   % Control
            data_mIL22{4:6,idx};...   % mIL22
            data_mIL22{7:9,idx};...   % LPS_10ng
            data_mIL22{10:12,idx};... % mIL22 + LPS_10ng
            data_mIL22{13:15,idx};... % LPS_100ng
            data_mIL22{16:18,idx};};  % mIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_mIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng mIL22';'10ng LPS';'1ng mIL22 + 10ng LPS';'100ng LPS';'1ng mIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('IL8 Response to mIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale),'.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/mIL22_boxplot_noVehicle_IL8','epsc')

% mIL22 MCP-1
idx = 11;
disp(data_mIL22.Properties.VariableNames(idx))
scale = mean(data_mIL22{1:3,idx});
plotData = {   % Control
            data_mIL22{4:6,idx};...   % mIL22
            data_mIL22{7:9,idx};...   % LPS_10ng
            data_mIL22{10:12,idx};... % mIL22 + LPS_10ng
            data_mIL22{13:15,idx};... % LPS_100ng
            data_mIL22{16:18,idx};};  % mIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_mIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng mIL22';'10ng LPS';'1ng mIL22 + 10ng LPS';'100ng LPS';'1ng mIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('MCP-1 Response to mIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/mIL22_boxplot_noVehicle_MCP1','epsc')

% mIL22 VEGF
idx = 5;
disp(data_mIL22.Properties.VariableNames(idx))
scale = mean(data_mIL22{1:3,idx});
plotData = {   % Control
            data_mIL22{4:6,idx};...   % mIL22
            data_mIL22{7:9,idx};...   % LPS_10ng
            data_mIL22{10:12,idx};... % mIL22 + LPS_10ng
            data_mIL22{13:15,idx};... % LPS_100ng
            data_mIL22{16:18,idx};};  % mIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_mIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng mIL22';'10ng LPS';'1ng mIL22 + 10ng LPS';'100ng LPS';'1ng mIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('VEGF Response to mIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/mIL22_boxplot_noVehicle_VEGF','epsc')

% mIL22 IP-10
idx = 6;
disp(data_mIL22.Properties.VariableNames(idx))
scale = mean(data_mIL22{1:3,idx});
plotData = {   % Control
            data_mIL22{4:6,idx};...   % mIL22
            data_mIL22{7:9,idx};...   % LPS_10ng
            data_mIL22{10:12,idx};... % mIL22 + LPS_10ng
            data_mIL22{13:15,idx};... % LPS_100ng
            data_mIL22{16:18,idx};};  % mIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_mIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng mIL22';'10ng LPS';'1ng mIL22 + 10ng LPS';'100ng LPS';'1ng mIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('IP-10 Response to mIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/mIL22_boxplot_noVehicle_IP-10','epsc')

% mIL22 PDGF-BB
idx = 17;
disp(data_mIL22.Properties.VariableNames(idx))
scale = mean(data_mIL22{1:3,idx});
plotData = {   % Control
            data_mIL22{4:6,idx};...   % mIL22
            data_mIL22{7:9,idx};...   % LPS_10ng
            data_mIL22{10:12,idx};... % mIL22 + LPS_10ng
            data_mIL22{13:15,idx};... % LPS_100ng
            data_mIL22{16:18,idx};};  % mIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_mIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng mIL22';'10ng LPS';'1ng mIL22 + 10ng LPS';'100ng LPS';'1ng mIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('PDGF-BB Response to mIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/mIL22_boxplot_noVehicle_PDGF-BB','epsc')

% mIL22 G-CSF
idx = 19;
disp(data_mIL22.Properties.VariableNames(idx))
scale = mean(data_mIL22{1:3,idx});
plotData = {   % Control
            data_mIL22{4:6,idx};...   % mIL22
            data_mIL22{7:9,idx};...   % LPS_10ng
            data_mIL22{10:12,idx};... % mIL22 + LPS_10ng
            data_mIL22{13:15,idx};... % LPS_100ng
            data_mIL22{16:18,idx};};  % mIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_mIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng mIL22';'10ng LPS';'1ng mIL22 + 10ng LPS';'100ng LPS';'1ng mIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('G-CSF Response to mIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/mIL22_boxplot_noVehicle_G-CSF','epsc')



% mIL22 IL1RA
idx = 21;
disp(data_mIL22.Properties.VariableNames(idx))
scale = mean(data_mIL22{1:3,idx});
plotData = {   % Control
            data_mIL22{4:6,idx};...   % mIL22
            data_mIL22{7:9,idx};...   % LPS_10ng
            data_mIL22{10:12,idx};... % mIL22 + LPS_10ng
            data_mIL22{13:15,idx};... % LPS_100ng
            data_mIL22{16:18,idx};};  % mIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_mIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng mIL22';'10ng LPS';'1ng mIL22 + 10ng LPS';'100ng LPS';'1ng mIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('IL1RA Response to mIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/mIL22_boxplot_noVehicle_IL1RA','epsc')

% mIL22 IL17
idx = 24;
disp(data_mIL22.Properties.VariableNames(idx))
scale = mean(data_mIL22{1:3,idx});
plotData = {   % Control
            data_mIL22{4:6,idx};...   % mIL22
            data_mIL22{7:9,idx};...   % LPS_10ng
            data_mIL22{10:12,idx};... % mIL22 + LPS_10ng
            data_mIL22{13:15,idx};... % LPS_100ng
            data_mIL22{16:18,idx};};  % mIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_mIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng mIL22';'10ng LPS';'1ng mIL22 + 10ng LPS';'100ng LPS';'1ng mIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('IL17 Response to mIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/mIL22_boxplot_noVehicle_IL17','epsc')
%----------------------------------------------------------------------------------------------------%
% Specific bar plots of cytokines
%{

geIL22 PLSDA Significant Cytokiens
- IL17, IP_10, VEGF, MIP_1alfa, MCP_1, IL12, IL_1ra

%}
%----------------------------------------------------------------------------------------------------%
% geIL22 IL17
idx = 20;
disp(data_geIL22.Properties.VariableNames(idx))
scale = mean(data_geIL22{1:3,idx});
plotData = {   % Control
            data_geIL22{4:6,idx};...   % geIL22
            data_geIL22{10:12,idx};... % LPS_10ng
            data_geIL22{7:9,idx};...    % geIL22 + LPS_10ng
            data_geIL22{16:18,idx};... % LPS_100ng
            data_geIL22{13:15,idx};};  % geIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_geIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng geIL22';'10ng LPS';'1ng geIL22 + 10ng LPS';'100ng LPS';'1ng geIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('IL17 Response to geIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/geIL22_boxplot_noVehicle_IL17','epsc')

% geIL22 IP10
idx = 12;
disp(data_geIL22.Properties.VariableNames(idx))
scale = mean(data_geIL22{1:3,idx});
plotData = {   % Control
            data_geIL22{4:6,idx};...   % geIL22
            data_geIL22{10:12,idx};... % LPS_10ng
            data_geIL22{7:9,idx};...    % geIL22 + LPS_10ng
            data_geIL22{16:18,idx};... % LPS_100ng
            data_geIL22{13:15,idx};};  % geIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_geIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng geIL22';'10ng LPS';'1ng geIL22 + 10ng LPS';'100ng LPS';'1ng geIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('IP10 Response to geIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/geIL22_boxplot_noVehicle_IP10','epsc')

% geIL22 VEGF
idx = 15;
disp(data_geIL22.Properties.VariableNames(idx))
scale = mean(data_geIL22{1:3,idx});
plotData = {   % Control
            data_geIL22{4:6,idx};...   % geIL22
            data_geIL22{10:12,idx};... % LPS_10ng
            data_geIL22{7:9,idx};...    % geIL22 + LPS_10ng
            data_geIL22{16:18,idx};... % LPS_100ng
            data_geIL22{13:15,idx};};  % geIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_geIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng geIL22';'10ng LPS';'1ng geIL22 + 10ng LPS';'100ng LPS';'1ng geIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('VEGF Response to geIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/geIL22_boxplot_noVehicle_VEGF','epsc')

% geIL22 MIP1A
idx = 23;
disp(data_geIL22.Properties.VariableNames(idx))
scale = mean(data_geIL22{1:3,idx});
plotData = {   % Control
            data_geIL22{4:6,idx};...   % geIL22
            data_geIL22{10:12,idx};... % LPS_10ng
            data_geIL22{7:9,idx};...    % geIL22 + LPS_10ng
            data_geIL22{16:18,idx};... % LPS_100ng
            data_geIL22{13:15,idx};};  % geIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_geIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng geIL22';'10ng LPS';'1ng geIL22 + 10ng LPS';'100ng LPS';'1ng geIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title('MIP1A Response to geIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/geIL22_boxplot_noVehicle_MIP1A','epsc')

% geIL22 MCP1
idx = 13;
disp(data_geIL22.Properties.VariableNames(idx))
scale = mean(data_geIL22{1:3,idx});
plotData = {   % Control
            data_geIL22{4:6,idx};...   % geIL22
            data_geIL22{10:12,idx};... % LPS_10ng
            data_geIL22{7:9,idx};...    % geIL22 + LPS_10ng
            data_geIL22{16:18,idx};... % LPS_100ng
            data_geIL22{13:15,idx};};  % geIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_geIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng geIL22';'10ng LPS';'1ng geIL22 + 10ng LPS';'100ng LPS';'1ng geIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title(' MCP1 Response to geIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/geIL22_boxplot_noVehicle_MCP1','epsc')

% geIL22 IL12
idx = 7;
disp(data_geIL22.Properties.VariableNames(idx))
scale = mean(data_geIL22{1:3,idx});
plotData = {   % Control
            data_geIL22{4:6,idx};...   % geIL22
            data_geIL22{10:12,idx};... % LPS_10ng
            data_geIL22{7:9,idx};...    % geIL22 + LPS_10ng
            data_geIL22{16:18,idx};... % LPS_100ng
            data_geIL22{13:15,idx};};  % geIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_geIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng geIL22';'10ng LPS';'1ng geIL22 + 10ng LPS';'100ng LPS';'1ng geIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title(' IL12 Response to geIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/geIL22_boxplot_noVehicle_IL12','epsc')


% geIL22 IL1RA
idx = 26;
disp(data_geIL22.Properties.VariableNames(idx))
scale = mean(data_geIL22{1:3,idx});
plotData = {   % Control
            data_geIL22{4:6,idx};...   % geIL22
            data_geIL22{10:12,idx};... % LPS_10ng
            data_geIL22{7:9,idx};...    % geIL22 + LPS_10ng
            data_geIL22{16:18,idx};... % LPS_100ng
            data_geIL22{13:15,idx};};  % geIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_geIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng geIL22';'10ng LPS';'1ng geIL22 + 10ng LPS';'100ng LPS';'1ng geIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title(' IL1RA Response to geIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/geIL22_boxplot_noVehicle_IL1RA','epsc')

% geIL22 IL8
idx = 10;
disp(data_geIL22.Properties.VariableNames(idx))
scale = mean(data_geIL22{1:3,idx});
plotData = {   % Control
            data_geIL22{4:6,idx};...   % geIL22
            data_geIL22{10:12,idx};... % LPS_10ng
            data_geIL22{7:9,idx};...    % geIL22 + LPS_10ng
            data_geIL22{16:18,idx};... % LPS_100ng
            data_geIL22{13:15,idx};};  % geIL22 + LPS_100ng

figure
h = boxplot(log2(cell2mat(plotData)./scale),data_geIL22.Condition(4:end)); % 
set(h, 'linewidth' ,1)
set(gca,'XTickLabel', {'1ng geIL22';'10ng LPS';'1ng geIL22 + 10ng LPS';'100ng LPS';'1ng geIL22 + 100ng LPS'})
hold on
xCenter = 1:numel(plotData); 
title(' IL8 Response to geIL22')
spread = 0.5; % 0=no spread; 0.5=random spread within box bounds (can be any value)
for i = 1:numel(plotData)
    plot(rand(size(plotData{i}))*spread -(spread/2) + xCenter(i), log2(plotData{i}./scale), '.k','linewidth', 3,'MarkerSize',40,'MarkerFaceColor','k')
end
saveas(gcf,'V3_Results/geIL22_boxplot_noVehicle_IL8','epsc')
%----------------------------------------------------------------------------------------------------%
% Functions
%----------------------------------------------------------------------------------------------------%

function [vipScores] = pls_vip(stats,XS,XL,YL,i)
% Calculate normalized PLS weights
W0 = bsxfun(@rdivide,stats.W(:,i),sqrt(sum(stats.W(:,i).^2,1)));
sumSq = sum(XS.^2,1).*sum(YL.^2,1);
% Calculate VIP scores for NCOMP components
vipScores = sqrt(size(XL,1) * sum(bsxfun(@times,sumSq,W0.^2),2) ./ sum(sumSq,2));
end
