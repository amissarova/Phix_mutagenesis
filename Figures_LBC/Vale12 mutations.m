%% load data
PDIR = '~/Develop/Phix_mutagenesis/' ; 
load([PDIR '/Data/T.mat']);

% remove unneeded rows
idx = ~T.IS_REF_ALLELE & T.NtSubstitute~='a'  & T.NtSubstitute~='c'   & T.NtSubstitute~='t'   & T.NtSubstitute~='g' ;
T = T(idx,:);

load Data/ANNO_SUBSTITUTION.mat
idx = ~ANNO_SUBSTITUTION.sameNt & ANNO_SUBSTITUTION.NtSubstitute~='a'  & ANNO_SUBSTITUTION.NtSubstitute~='c'   & ANNO_SUBSTITUTION.NtSubstitute~='t'   & ANNO_SUBSTITUTION.NtSubstitute~='g' ;
ANNO_SUBSTITUTION = ANNO_SUBSTITUTION( idx , : );
%%
fh = figure('units','centimeters','position',[5 5 8 8]);
hold on ;
position = 4972 ; 
xl = linspace(0,5,1e3);
for L = 'ACTG'
    dataA = T.sub_prentnt_z_mode(   T.NtSubstitute == L)  ; 
    data = T.sub_prentnt_z_mode( T.Position == position & T.NtSubstitute == L) ; 
    y = NaN(numel(xl),2);
    for I = 1:numel(xl)
        y(I,1) = 100*nanmean(data>xl(I));
        y(I,2) = 100*nanmean(dataA>xl(I));
    end
    if ~all(isnan(y(:,1)))
        plot(xl,y(:,1),'-','DisplayName',L,'LineWidth',2);
    end
    plot(xl,y(:,2),'-','DisplayName',['all ' L],'LineWidth',2,'Color',[.7 .7 .7]);

end
legend('location','best')
title(['Position = ' num2str(position) ] );
ylabel('% of expts w/z > X')
xlabel('Z threshold')
axis tight; 
ylim([0 50])
set(gca,'xtick',-10:10)
grid on ;