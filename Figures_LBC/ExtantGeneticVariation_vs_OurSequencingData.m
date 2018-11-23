%% how do our observed allele frequencies compare to those found in natural populations? 
% 
% Questions to ask: 
%
% Do alleles found in nature have higher Z scores? 
% Do alleles NOT found in nature have negative Z scores? 
% 
% Within a set alleles groups by Zscore, how many around found in nature?
%    are more found in nature in groupings at high Z ? 

%% load sequencing data
PDIR = '~/Develop/Phix_mutagenesis/'  ; 
load([ PDIR 'Data/T.mat' ]);

%% load natural population data
NatPop = AlignMultipleFASTAs_ProduceTabOfVariation( ) ; 

%% Take mean & median z-score across all experiments
% and join the two datasets
T.zGT3 = T.sub_prentnt_z_mode > 3 ; 
T.zGT2p5 = T.sub_prentnt_z_mode > 2.5 ; 
T.zGT2 = T.sub_prentnt_z_mode > 2 ; 
T.zGT1 = T.sub_prentnt_z_mode > 1 ; 
T.zGT1p5 = T.sub_prentnt_z_mode > 1.5 ; 
T.zLT1   = T.sub_prentnt_z_mode < -1 ; 
T.zLT1p5   = T.sub_prentnt_z_mode < -1.5 ; 
T.zLT2   = T.sub_prentnt_z_mode < -2 ; 
T.zLT2p5   = T.sub_prentnt_z_mode < -2.5 ; 

vn = T.Properties.VariableNames ; 

G = grpstats( T , {'Position' 'Nt' 'NtSubstitute'} , {'median' 'mean' 'std' 'sem'} , 'DataVars' , vn(regexpcmp(vn,'z')) );

G = innerjoin( G , NatPop , 'Key', {'Position' 'Nt' 'NtSubstitute'});
G = G( G.Nt ~= G.NtSubstitute , :);

%% join T & NatPop
G2 = outerjoin( T(T.Nt ~= T.NtSubstitute & T.NtSubstitute~='a'& T.NtSubstitute~='c'& T.NtSubstitute~='t'& T.NtSubstitute~='g',:) , NatPop(NatPop.NtSubstitute ~= 'N' & NatPop.NtSubstitute ~= '-' & NatPop.NtSubstitute ~= NatPop.Nt , :) , 'Key', {'Position' 'Nt' 'NtSubstitute'});

%% ecdfs 
fh = figure('units','centimeters','position',[5 5 18 18]);
subplot(2,2,1)
hold on ; 
[f,x] = ecdf( G.median_sub_prentnt_z_mode( G.AlleleN==0));
plot(x,f,'-k','LineWidth',2,'DisplayName','never occurs')
[f,x] = ecdf( G.median_sub_prentnt_z_mode( G.AlleleN>0));
plot(x,f,'-','LineWidth',2,'DisplayName','occurs >= 1')
xlim([-2 5])
legend('location','nw')
xlabel('median Z score across experiments')
ylabel('Fraction of positions in genome')
title('median Z')

subplot(2,2,2)
hold on ; 
[f,x] = ecdf( G.mean_sub_prentnt_z_mode( G.AlleleN==0));
plot(x,f,'-k','LineWidth',2,'DisplayName','never occurs')
[f,x] = ecdf( G.mean_sub_prentnt_z_mode( G.AlleleN>0));
plot(x,f,'-','LineWidth',2,'DisplayName','occurs >= 1')
xlim([-2 5])
legend('location','nw')
xlabel('mean Z score across experiments')
ylabel('Fraction of positions in genome')
title('mean Z')


subplot(2,2,3)
hold on ; 
[f,x] = ecdf( G.mean_zGT2( G.AlleleN==0));
plot(x,f,'-k','LineWidth',2,'DisplayName','never occurs')
[f,x] = ecdf( G.mean_zGT2( G.AlleleN>0));
plot(x,f,'-','LineWidth',2,'DisplayName','occurs >= 1')
xlim([0 1])
legend('location','se')
xlabel('Fraction of expts Z > 2')
ylabel('Fraction of positions in genome')
title('Z > 2')


subplot(2,2,4)
hold on ; 
[f,x] = ecdf( G.mean_zGT3( G.AlleleN==0));
plot(x,f,'-k','LineWidth',2,'DisplayName','never occurs')
[f,x] = ecdf( G.mean_zGT3( G.AlleleN>0));
plot(x,f,'-','LineWidth',2,'DisplayName','occurs >= 1')
xlim([0 1])
legend('location','se')
xlabel('Fraction of expts Z > 3')
ylabel('Fraction of positions in genome')
title('Z > 3')

%% histograms 
xl = linspace(-2,5,50);

fh = figure('units','centimeters','position',[5 5 18 18]);
subplot(2,2,1)
y = G.median_sub_prentnt_z_mode; y(y>5)=5 ; y(y<-2)=-2 ; 
hold on ; 
histogram( y( G.AlleleN==0) , xl , 'Normalization','Probability');
histogram( y( G.AlleleN>0) , xl , 'Normalization','Probability');
xlim([-2 5.2])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('median Z score across experiments')
ylabel('Fraction of positions in genome')
title('median Z')


subplot(2,2,2)
hold on ; 
y = G.mean_sub_prentnt_z_mode; y(y>5)=5 ; y(y<-2)=-2 ; 
histogram( y( G.AlleleN==0) , 50 , 'Normalization','Probability');
histogram( y( G.AlleleN>0) , 50 , 'Normalization','Probability');
xlim([-2 5.2])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('mean Z score across experiments')
ylabel('Fraction of positions in genome')
title('mean Z')

xl = linspace(0,1,50);

subplot(2,2,3)
hold on ; 
histogram( G.mean_zGT2( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zGT2( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z > 2')
ylabel('Fraction of positions in genome')
title('Z > 2')

subplot(2,2,4)
hold on ; 
histogram( G.mean_zGT3( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zGT3( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z > 3')
ylabel('Fraction of positions in genome')
title('Z > 3')



%% %% %% %% %% histograms continuous Z > X 1-2.5

xl = linspace(0,1,50);
fh = figure('units','centimeters','position',[5 5 18 18]);

subplot(2,2,1)
hold on ; 
histogram( G.mean_zGT1( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zGT1( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z > 1')
ylabel('Fraction of positions in genome')
title('Z > 1')

subplot(2,2,2)
hold on ; 
histogram( G.mean_zGT1p5( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zGT1p5( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z > 1.5')
ylabel('Fraction of positions in genome')
title('Z > 1.5')

subplot(2,2,3)
hold on ; 
histogram( G.mean_zGT2( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zGT2( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z > 2')
ylabel('Fraction of positions in genome')
title('Z > 2')

subplot(2,2,4)
hold on ; 
histogram( G.mean_zGT2p5( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zGT2p5( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z > 2.5')
ylabel('Fraction of positions in genome')
title('Z > 2.5')


%%
fh = figure('units','centimeters','position',[5 5 18 18]);

subplot(2,2,1)
hold on ; 
histogram( G.mean_zLT1( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zLT1( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z < -1')
ylabel('Fraction of positions in genome')
title('Z < -1')

subplot(2,2,2)
hold on ; 
histogram( G.mean_zLT1p5( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zLT1p5( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z < -1.5')
ylabel('Fraction of positions in genome')
title('Z < -1.5')

subplot(2,2,3)
hold on ; 
histogram( G.mean_zLT2( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zLT2( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z < -2')
ylabel('Fraction of positions in genome')
title('Z < -2')
ylim([0 0.05])

subplot(2,2,4)
hold on ; 
histogram( G.mean_zLT2p5( G.AlleleN==0),xl , 'Normalization','Probability');
histogram( G.mean_zLT2p5( G.AlleleN>0),xl , 'Normalization','Probability');
xlim([0 1])
legend({'never occurs' 'occurs >= 1'}, 'location','nw')
xlabel('Fraction of expts Z < -2.5')
ylabel('Fraction of positions in genome')
title('Z < -2.5')
ylim([0 0.05])
%%
xl = linspace( 0,10 , 25 );
y = NaN( numel(xl) , 1);
for I = 1:numel(xl)
    idx_data_gt_thresh = G2.sub_prentnt_z_mode > xl(I)-0.25 &  G2.sub_prentnt_z_mode < xl(I)+0.25; 
    idx_found_in_nature_0 = G2.AlleleN>0  ;
    y(I,1) = sum(idx_data_gt_thresh & idx_found_in_nature_0)  / sum(idx_data_gt_thresh) ; 
end
fh = figure('units','centimeters','position',[5 5 10 10]);
plot(xl,100*y,'ok','MarkerFaceColor',[.7 .7 .7])
xlabel('Z score')
set(gca,'xtick',-10:10)
ylabel('% of alleles with this Z score that are found in nature')
%%  hmm. bug?
%id = strcat( string(G2.PreNt) , string( G2.Nt_left) , string(G2.NtSubstitute_left)) ; 
xl = linspace( 0,10 , 50 );

id = strcat( string( G2.Nt_left) , string(G2.NtSubstitute_left)) ; 

uid = unique(id) ;

y2 = NaN( numel(xl) , numel(uid));
for idI = 1:numel(uid)
    idx = strcmp(id,uid{idI});
    for I = 1:numel(xl)
        idx_data_gt_thresh = G2.sub_prentnt_z_mode > xl(I)-0.5 &  G2.sub_prentnt_z_mode < xl(I)+0.5; 
        idx_found_in_nature_0 = G2.AlleleN>0  ;
        y2(I,idI) = sum(idx_data_gt_thresh & idx_found_in_nature_0 & idx)  / sum(idx_data_gt_thresh & idx) ; 
    end
end

figure; 
hold on ; 
set(gca,'ColorOrder',parula(numel(uid))) ; 
plot(xl,(y2'),'.-'); 
legend(uid);
xlabel('Z score')
set(gca,'xtick',-10:10)
ylabel('% of alleles with this Z score that are found in nature')
%%
