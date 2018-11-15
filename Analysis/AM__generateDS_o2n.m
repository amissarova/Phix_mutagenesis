%% 1. load template DS_MutationsInfo_Illumina

cd ~/Develop/Phix_mutagenesis/
addpath(genpath('~/Develop/matlab'));
load('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat'); 
D = DS;

%% 2. add Frequencies and Coverage for o2n data
samples = dir('~/Develop/Phix_mutagenesis/o2n_Data'); 
cd ~/Develop/Phix_mutagenesis/o2n_Data
for I = 1:length(samples)
    if findstr(samples(I).name , 'adp')
        fprintf('Processed %.2f.\n' , I/length(samples));
        T = dataset('file' , samples(I).name);
        T.Properties.VarNames{2} = 'Position';
        D.CovW = NaN(length(D) , 1);
        D.FreqW = NaN(length(D) , 1);
        D.CovC = NaN(length(D) , 1);
        D.FreqC = NaN(length(D) , 1);
        D.Freq = NaN(length(D) , 1);
        for J = 1:length(D)
            idx = find(D.Position(J) == T.Position);
            names_vars_T = get(T , 'VarNames');
            idx_substitute_W = find(strcmp(D.NtSubstitute{J} , names_vars_T));
            idx_substitute_C = find(strcmp(lower(D.NtSubstitute{J}) , names_vars_T));
            if T.covW(idx) > 0
                D.CovW(J) = double(T(idx , idx_substitute_W));
                D.FreqW(J) = double(T(idx , idx_substitute_W))/T.covW(idx);
            end
            if T.covC(idx) > 0
                D.CovC(J) = double(T(idx , idx_substitute_C));
                D.FreqC(J) = double(T(idx , idx_substitute_C))/T.covC(idx);
            end
            if T.cov(idx) > 0
                D.Freq(J) = (double(T(idx , idx_substitute_C)) + double(T(idx , idx_substitute_W)))/T.cov(idx);
            end
        end
        names_vars_D = get(D , 'VarNames');
        
        prefix_sample = strrep(samples(I).name , '-trim.srt_RG.bam.mpileup.tab' , '');
        prefix_sample = strrep(prefix_sample , 'adp-' , '');
        idx = find(strcmp(names_vars_D , 'CovW'));
        D.Properties.VarNames{idx} = strcat('CovW_' , prefix_sample);
        idx = find(strcmp(names_vars_D , 'FreqW'));
        D.Properties.VarNames{idx} = strcat('FreqW_' , prefix_sample);
        idx = find(strcmp(names_vars_D , 'CovC'));
        D.Properties.VarNames{idx} = strcat('CovC_' , prefix_sample);
        idx = find(strcmp(names_vars_D , 'FreqC'));
        D.Properties.VarNames{idx} = strcat('FreqC_' , prefix_sample);
        idx = find(strcmp(names_vars_D , 'Freq'));
        D.Properties.VarNames{idx} = strcat('Freq_' , prefix_sample);
    end
end
cd ..

%
names_vars_D = get(D , 'VarNames');

D.vecFreqW_o2n = cell(length(D) , 1);
D.vecCovW_o2n = cell(length(D) , 1);
D.vecFreqC_o2n = cell(length(D) , 1);
D.vecCovC_o2n = cell(length(D) , 1);
D.vecFreq_o2n = cell(length(D) , 1);

idx = [];
for I = 1:length(names_vars_D)
    if findstr('FreqW_' , names_vars_D{I})
        idx = [idx I];
    end
end
for I = 1:length(D)
    vecFrequency = NaN(length(idx) , 1);
    for J = 1:length(idx)
        vecFrequency(J) = D(I , idx(J));
    end
    D.vecFreqW_o2n{I} = vecFrequency;
end
idx = [];
for I = 1:length(names_vars_D)
    if findstr('FreqC_' , names_vars_D{I})
        idx = [idx I];
    end
end
for I = 1:length(D)
    vecFrequency = NaN(length(idx) , 1);
    for J = 1:length(idx)
        vecFrequency(J) = D(I , idx(J));
    end
    D.vecFreqC_o2n{I} = vecFrequency;
end
idx = [];
for I = 1:length(names_vars_D)
    if findstr('CovW_' , names_vars_D{I})
        idx = [idx I];
    end
end
for I = 1:length(D)
    vecFrequency = NaN(length(idx) , 1);
    for J = 1:length(idx)
        vecFrequency(J) = D(I , idx(J));
    end
    D.vecCovW_o2n{I} = vecFrequency;
end
idx = [];
for I = 1:length(names_vars_D)
    if findstr('CovC_' , names_vars_D{I})
        idx = [idx I];
    end
end
for I = 1:length(D)
    vecFrequency = NaN(length(idx) , 1);
    for J = 1:length(idx)
        vecFrequency(J) = D(I , idx(J));
    end
    D.vecCovC_o2n{I} = vecFrequency;
end
idx = [];
for I = 1:length(names_vars_D)
    if findstr('Freq_' , names_vars_D{I})
        idx = [idx I];
    end
end
for I = 1:length(D)
    vecFrequency = NaN(length(idx) , 1);
    for J = 1:length(idx)
        vecFrequency(J) = D(I , idx(J));
    end
    D.vecFreq_o2n{I} = vecFrequency;
end

%
D.medianFreqW_o2n = NaN(length(D) , 1);
D.medianFreqC_o2n = NaN(length(D) , 1);
D.medianFreq_o2n = NaN(length(D) , 1);
for I = 1:length(D)
	D.medianFreqW_o2n(I) = nanmedian(D.vecFreqW_o2n{I});
    D.medianFreqC_o2n(I) = nanmedian(D.vecFreqC_o2n{I});
    D.medianFreq_o2n(I) = nanmedian(D.vecFreq_o2n{I});
end
D.medianCovW_o2n = NaN(length(D) , 1);
D.medianCovC_o2n = NaN(length(D) , 1);
for I = 1:length(D)
	D.medianCovW_o2n(I) = nanmedian(D.vecCovW_o2n{I})./nanmedian(D.vecFreqW_o2n{I});
    D.medianCovC_o2n(I) = nanmedian(D.vecCovC_o2n{I})./nanmedian(D.vecFreqC_o2n{I});
end
save('~/Develop/Phix_Mutagenesis/Data/DS__o2n.mat' , 'D'); 



    
    
    
