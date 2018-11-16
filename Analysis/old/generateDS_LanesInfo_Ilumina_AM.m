%% 1. Pull the template DS_mutations_info w/ info regarding every possible substitution

cd ~/Develop/Phix_mutagenesis/
addpath(genpath('~/Develop/matlab'));
load('~/Develop/Phix_Mutagenesis/Data/DS_MutationsInfo_Ilumina.mat'); 
D = DS;
D = D(: , {'Position', 'Nt', 'NtSubstitute' , 'PreNtNtNtSub'});
% 2. DS w lane info for the frist DS
samples = dir('~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane1'); 
cd ~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane1
for I = 1:length(samples)
    if findstr(samples(I).name , 'StrandsSep')
        fprintf('Processed %.2f.\n' , I/length(samples));
        % read each dataset
        T = dataset('file' , samples(I).name);
        T.Properties.VarNames{2} = 'Position';

        % Cov -- FWD and REV
        D.CovW_temp = NaN(length(D) , 1);
        D.CovC_temp = NaN(length(D) , 1);

        % Freq -- FWD, REV and TOTAL
        D.FreqW_temp = NaN(length(D) , 1);
        D.FreqC_temp = NaN(length(D) , 1);
        D.Freq_temp = NaN(length(D) , 1);

        for J = 1:length(D)
            idx = find(D.Position(J) == T.Position);
            names_vars_T = get(T , 'VarNames');
            idx_substitute_W = find(strcmp(D.NtSubstitute{J} , names_vars_T));
            idx_substitute_C = find(strcmp(lower(D.NtSubstitute{J}) , names_vars_T));
            if T.covW(idx) > 0
                D.CovW_temp(J) = double(T(idx , idx_substitute_W));
                D.FreqW_temp(J) = double(T(idx , idx_substitute_W))/T.covW(idx);
            end
            if T.covC(idx) > 0
                D.CovC_temp(J) = double(T(idx , idx_substitute_C));
                D.FreqC_temp(J) = double(T(idx , idx_substitute_C))/T.covC(idx);
            end
            if T.cov(idx) > 0
                D.Freq_temp(J) = (double(T(idx , idx_substitute_C)) + double(T(idx , idx_substitute_W)))/T.cov(idx);
            end
        end
        names_vars_D = get(D , 'VarNames');
        
        prefix_sample = strrep(samples(I).name , '.pileupCountAllAlleles_StrandsSep.tab' , '');
        prefix_sample = strrep(prefix_sample , 'NoIndex_L00' , '');
        prefix_sample = strrep(prefix_sample , '_7001450_0407_ACC8GKANXX_' , '');
        prefix_sample = strrep(prefix_sample , '_7001450_0415_ACC843ANXX_' , '');
        prefix_sample = strrep(prefix_sample , 'lane' , '');
        idx = find(strcmp(names_vars_D , 'CovW_temp'));
        D.Properties.VarNames{idx} = strcat('CovW_' , prefix_sample);
        idx = find(strcmp(names_vars_D , 'FreqW_temp'));
        D.Properties.VarNames{idx} = strcat('FreqW_' , prefix_sample);
        idx = find(strcmp(names_vars_D , 'CovC_temp'));
        D.Properties.VarNames{idx} = strcat('CovC_' , prefix_sample);
        idx = find(strcmp(names_vars_D , 'FreqC_temp'));
        D.Properties.VarNames{idx} = strcat('FreqC_' , prefix_sample);
        idx = find(strcmp(names_vars_D , 'Freq_temp'));
        D.Properties.VarNames{idx} = strcat('Freq_' , prefix_sample);
    end
end
cd ..

% 3. combine all lanes together
names_vars_D = get(D , 'VarNames');

D.vecFreqW = cell(length(D) , 1);
D.vecCovW = cell(length(D) , 1);
D.vecFreqC = cell(length(D) , 1);
D.vecCovC = cell(length(D) , 1);
D.vecFreq = cell(length(D) , 1);

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
    D.vecFreqW{I} = vecFrequency;
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
    D.vecFreqC{I} = vecFrequency;
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
    D.vecCovW{I} = vecFrequency;
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
    D.vecCovC{I} = vecFrequency;
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
    D.vecFreq{I} = vecFrequency;
end
%
%D = D(: , {'Position', 'Nt', 'NtSubstitute' , ...
%    'vecFreqW','vecCovW','vecFreqC','vecCovC','vecFreq'});

% 4. add median across all lanes
D.medianFreqW = cell(length(D) , 1);
D.medianFreqC = cell(length(D) , 1);
D.medianFreq = cell(length(D) , 1);
D.medianCovW = cell(length(D) , 1);
D.medianCovC = cell(length(D) , 1);
for I = 1:length(D)
	D.medianFreqW{I} = [nanmedian(D.vecFreqW{I})];
    D.medianFreqC{I} = [nanmedian(D.vecFreqC{I})];
    D.medianFreq{I} = [nanmedian(D.vecFreq{I})];
	D.medianCovW{I} = [nanmedian(D.vecCovW{I})];
    D.medianCovC{I} = [nanmedian(D.vecCovC{I})];
end
D = D(: , {'Position', 'Nt', 'NtSubstitute' , 'PreNtNtNtSub' , ...
    'medianFreqW','medianCovW','medianFreqC','medianCovC','medianFreq'});

% 5. add experiment ID
D.experimentID = cell(length(D) , 1);
ID = cell(1); ID{1} = '180316_7001450_0407_ACC843ANXX__lane1';
for I = 1:length(D)
    D.experimentID{I} = ID;
end

% 6. add mode, std and z-score
%load('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat');
D.m_mode = cell(length(D) , 1);
D.s_mode = cell(length(D) , 1);
D.z_mode = cell(length(D) , 1);

unq_id = unique(D.PreNtNtNtSub); unq_id = setdiff(unq_id , {'same'}); unq_id = sort(unq_id);
for J = 1:length(unq_id)
	
    data = NaN(length(D) , 1);
    for I = 1:length(D)
        temp_data = D.medianFreqW{I};
        data(I) = temp_data(end);
    end
    idx = find(strcmp( D.PreNtNtNtSub , unq_id{J}) );
	T = D(idx , :);
    
    data_T = NaN(length(T) , 1);
    for I = 1:length(T)
        temp_data = T.medianFreqW{I};
        data_T(I) = temp_data(end);
    end

    % m and s
    data_T1 = data_T(~isnan(data_T)); data_T1 = data_T1(data_T1 < Inf);
	m_mode = modefit(data_T1 , 0 , [-0.000001:0.000001:0.005]);
    temp_data = data_T1(data_T1 <= m_mode);
	data = temp_data;
	for I = 1:length(temp_data)
        data = [data ; m_mode + m_mode - temp_data(I)];
    end    
	s_mode = nanstd(data);
    
    % 
	for Z = 1:length(idx)
        D.m_mode{idx(Z)} = [m_mode];
        D.s_mode{idx(Z)} = [s_mode];
        D.z_mode{idx(Z)} = [(data_T(Z) - m_mode)/s_mode];
    end
end

% 7. join w original DS
D = join(DS , D , 'Keys' , {'Position' , 'Nt' , 'NtSubstitute' , 'PreNtNtNtSub'} , 'Type' , 'Left' , 'MergeKeys' , true);
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
%%
% 8. add other Datasets
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane2');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
% BioRep-0316
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane3');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane4');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane5');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane6');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane7');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180316_7001450_0407_ACC843ANXX__lane8');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 

D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180402_7001450_0415_ACC843ANXX__1');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180402_7001450_0415_ACC843ANXX__2');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180402_7001450_0415_ACC843ANXX__3');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180402_7001450_0415_ACC843ANXX__4');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180402_7001450_0415_ACC843ANXX__5');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180402_7001450_0415_ACC843ANXX__6');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 
D = addPhixMutationIlluminaDataset(D , '~/Develop/Phix_mutagenesis/Illumina_Data/Data_StrandsSep/180402_7001450_0415_ACC843ANXX__8');
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 


%% OLD, but keep in case

%%
load('~/Develop/Phix_Mutagenesis/Data/DS_ConservationMtrxProteinA.mat');
load('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat');
%
% if position in A ORF
D.InA_bool = zeros(length(D) , 1);
% nonsense in A/nonstop in A/synonymous in A/missense-not
% conserved/missense conserved
D.InA_class = cell(length(D) , 1);
% for missense-conserved -- max score./1000
D.InA_value = zeros(length(D) , 1);
for I = 1:length(D)
    if ~strcmp(D.Nt{I} , D.NtSubstitute{I})
        genes = D.Genes{I};
        if ~isempty(find(strcmp(genes , 'A')))
            D.InA_bool(I) = 1;
            idx = find(strcmp(genes , 'A'));
            AA = D.AAs{I}; AA = AA(idx);
            AASubstitute = D.AASubstitutes{I}; AASubstitute = AASubstitute(idx);
            PositionInProtein = D.PositionInProteins{I}; PositionInProtein = PositionInProtein(idx);
            if strcmp(AASubstitute , 'STOP') & ~strcmp(AA , 'STOP')
                D.InA_class{I} = 'Nonsense';
            elseif strcmp(AA , 'STOP') & ~strcmp(AASubstitute , 'STOP')
                D.InA_class{I} = 'Nonstop';
            elseif strcmp(AA , AASubstitute)
                D.InA_class{I} = 'Synonymous';
            else
                idx = find(PositionInProtein == ConservationMtrx.PositionAA & ...
                    strcmp(AA , ConservationMtrx.AA) & strcmp(AASubstitute , ConservationMtrx.AASubstitute));
                if ~isempty(idx)
                    D.InA_class{I} = 'Conserved';
                    D.InA_value(I) = ConservationMtrx.max_Score(idx)./1000;
                else 
                    idx = find(PositionInProtein == ConservationMtrx.PositionAA & ...
                    strcmp(AA , ConservationMtrx.AA) & strcmp('-' , ConservationMtrx.AASubstitute));
                    if ~isempty(idx)
                        D.InA_class{I} = 'ConservedInDel';
                        D.InA_value(I) = ConservationMtrx.max_Score(idx)./1000;
                    else
                        D.InA_class{I} = 'NotConserved';
                    end

                end
            end
        end
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 

%%
load('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat');
D.m1_0316 = NaN(length(D) , 1);
D.m2_0316 = NaN(length(D) , 1);
D.m3_0316 = NaN(length(D) , 1);

D.m1_0402 = NaN(length(D) , 1);
D.m2_0402 = NaN(length(D) , 1);
D.m3_0402 = NaN(length(D) , 1);

D.s1_0316 = NaN(length(D) , 1);
D.s2_0316 = NaN(length(D) , 1);
D.s3_0316 = NaN(length(D) , 1);

D.s1_0402 = NaN(length(D) , 1);
D.s2_0402 = NaN(length(D) , 1);
D.s3_0402 = NaN(length(D) , 1);

D.z1_0316 = NaN(length(D) , 1);
D.z2_0316 = NaN(length(D) , 1);
D.z3_0316 = NaN(length(D) , 1);

D.z1_0402 = NaN(length(D) , 1);
D.z2_0402 = NaN(length(D) , 1);
D.z3_0402 = NaN(length(D) , 1);

D.LLM1_0316 = zeros(length(D) , 1);
D.LLM2_0316 = zeros(length(D) , 1);
D.LLM3_0316 = zeros(length(D) , 1);

D.LLM1_0402 = zeros(length(D) , 1);
D.LLM2_0402 = zeros(length(D) , 1);
D.LLM3_0402 = zeros(length(D) , 1);

unq_id = unique(D.PreNtNtNtSub); unq_id = setdiff(unq_id , {'same'}); unq_id = sort(unq_id);
for J = 1:length(unq_id)
	% 03-16
    idx = find(strcmp( D.PreNtNtNtSub , unq_id{J}) & ...
        D.medianFreqW_180316 < Inf);
	T = D(idx , :);
    
    % m1 and s1
	m1 = modefit(T.medianFreqW_180316  , 0 , [-0.000001:0.000001:0.005]);
    temp_data = T.medianFreqW_180316(T.medianFreqW_180316 <= m1);
	data = temp_data;
	for I = 1:length(temp_data)
        data = [data ; m1 + m1 - temp_data(I)];
    end    
	s1 = nanstd(data);
    
    % m2 and s2
    m2 = nanmedian(T.medianFreqW_180316);
    temp_data = T.medianFreqW_180316(T.medianFreqW_180316 <= m2);
	data = temp_data;
	for I = 1:length(temp_data)
        data = [data ; m2 + m2 - temp_data(I)];
    end    
	s2 = nanstd(data);
    
    % m3 and s3
	idx_A = find(T.InA_bool == 1);
    m3 = nanmedian(T.medianFreqW_180316(idx_A));
    temp_data = T.medianFreqW_180316(idx_A);
    temp_data = temp_data(temp_data <= m3);
	data = temp_data;
	for I = 1:length(temp_data)
        data = [data ; m3 + m3 - temp_data(I)];
    end    
	s3 = nanstd(data);
    
	for Z = 1:length(idx)
        D.m1_0316(idx(Z)) = m1;
        D.s1_0316(idx(Z)) = s1;
        D.z1_0316(idx(Z)) = (D.medianFreqW_180316(idx(Z)) - m1)/s1;
        if D.z1_0316(idx(Z)) > 3
            D.LLM1_0316(idx(Z)) = 1;
        end
        
        D.m2_0316(idx(Z)) = m2;
        D.s2_0316(idx(Z)) = s2;
        D.z2_0316(idx(Z)) = (D.medianFreqW_180316(idx(Z)) - m2)/s2;
        if D.z2_0316(idx(Z)) > 3
            D.LLM2_0316(idx(Z)) = 1;
        end
        
        D.m3_0316(idx(Z)) = m3;
        D.s3_0316(idx(Z)) = s3;
        D.z3_0316(idx(Z)) = (D.medianFreqW_180316(idx(Z)) - m3)/s3;
        if D.z3_0316(idx(Z)) > 3
            D.LLM3_0316(idx(Z)) = 1;
        end
    end
    % 04-02
    idx = find(strcmp( D.PreNtNtNtSub , unq_id{J}) & ...
        D.medianFreqW_180402 < Inf);
	T = D(idx , :);
    
    % m1 and s1
	m1 = modefit(T.medianFreqW_180402  , 0 , [-0.000001:0.000001:0.005]);
    temp_data = T.medianFreqW_180402(T.medianFreqW_180402 <= m1);
	data = temp_data;
	for I = 1:length(temp_data)
        data = [data ; m1 + m1 - temp_data(I)];
    end    
	s1 = nanstd(data);
    
    % m2 and s2
    m2 = nanmedian(T.medianFreqW_180402);
    temp_data = T.medianFreqW_180402(T.medianFreqW_180402 <= m2);
	data = temp_data;
	for I = 1:length(temp_data)
        data = [data ; m2 + m2 - temp_data(I)];
    end    
	s2 = nanstd(data);
    
    % m3 and s3
    idx_A = find(T.InA_bool == 1);
    m3 = nanmedian(T.medianFreqW_180402(idx_A));
    temp_data = T.medianFreqW_180402(idx_A);
    temp_data = temp_data(temp_data <= m3);
	data = temp_data;
	for I = 1:length(temp_data)
        data = [data ; m3 + m3 - temp_data(I)];
    end    
	s3 = nanstd(data);
    
	for Z = 1:length(idx)
        D.m1_0402(idx(Z)) = m1;
        D.s1_0402(idx(Z)) = s1;
        D.z1_0402(idx(Z)) = (D.medianFreqW_180402(idx(Z)) - m1)/s1;
        if D.z1_0402(idx(Z)) > 3
            D.LLM1_0402(idx(Z)) = 1;
        end
        
        D.m2_0402(idx(Z)) = m2;
        D.s2_0402(idx(Z)) = s2;
        D.z2_0402(idx(Z)) = (D.medianFreqW_180402(idx(Z)) - m2)/s2;
        if D.z2_0402(idx(Z)) > 3
            D.LLM2_0402(idx(Z)) = 1;
        end
        
        D.m3_0402(idx(Z)) = m3;
        D.s3_0402(idx(Z)) = s3;
        D.z3_0402(idx(Z)) = (D.medianFreqW_180402(idx(Z)) - m3)/s3;
        if D.z3_0402(idx(Z)) > 3
            D.LLM3_0402(idx(Z)) = 1;
        end
        
    end
end
save('~/Develop/Phix_Mutagenesis/Data/DS_LanesInfo_Ilumina.mat' , 'D'); 















    
    
    
