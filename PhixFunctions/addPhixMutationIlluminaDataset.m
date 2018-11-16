function D = addPhixMutationIlluminaDataset(DS , dataDir)
%%
% 1. pull the template
D = DS;
D = D(: , {'Position', 'Nt', 'NtSubstitute' , 'PreNtNtNtSub'});

% 2. DS w lane info for the frist DS
samples = dir(dataDir); 
cd(dataDir)
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
ID = cell(1); 
temp_ID = regexp( dataDir , '/' ,'split'); temp_ID = temp_ID{end};
temp_ID = strrep(temp_ID , '.tab' , '');
ID{1} = temp_ID;
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
    data_T1 = data_T(~isnan(data_T)); data_T1 = data_T1(data_T1 < Inf);
    % m and s
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

% 7. merge DSs
var_names_DS = get(DS , 'VarNames');
if isempty(find(strcmp(var_names_DS , 'medianFreqW')))
    DS.medianFreqW = cell(length(DS) , 1);
    DS.medianFreqC = cell(length(DS) , 1);
    DS.medianFreq = cell(length(DS) , 1);
    DS.medianCovW = cell(length(DS) , 1);
    DS.medianCovC = cell(length(DS) , 1);
    DS.m_mode = cell(length(DS) , 1);
    DS.s_mode = cell(length(DS) , 1);
    DS.z_mode = cell(length(DS) , 1);
    DS.experimentID = cell(length(DS) , 1);
    for I = 1:length(DS) 
        DS.medianFreqW{I} = [D.medianFreqW{I}]; 
        DS.medianFreqC{I} = [D.medianFreqC{I}]; 
        DS.medianFreq{I} = [D.medianFreq{I}]; 
        DS.medianCovW{I} = [D.medianCovW{I}]; 
        DS.medianCovC{I} = [D.medianCovC{I}]; 
        DS.m_mode{I} = [D.m_mode{I}]; 
        DS.s_mode{I} = [D.s_mode{I}];
        DS.z_mode{I} = [D.z_mode{I}];
        DS.experimentID{I} = D.experimentID{I};
    end
    
else
    for I = 1:length(DS)
        DS.medianFreqW{I} = [DS.medianFreqW{I} D.medianFreqW{I}]; 
        DS.medianFreqC{I} = [DS.medianFreqC{I} D.medianFreqC{I}]; 
        DS.medianFreq{I} = [DS.medianFreq{I} D.medianFreq{I}]; 
        DS.medianCovW{I} = [DS.medianCovW{I} D.medianCovW{I}]; 
        DS.medianCovC{I} = [DS.medianCovC{I} D.medianCovC{I}]; 
        DS.m_mode{I} = [DS.m_mode{I} D.m_mode{I}]; 
        DS.s_mode{I} = [DS.s_mode{I} D.s_mode{I}];
        DS.z_mode{I} = [DS.z_mode{I} D.z_mode{I}];
        experimentID = DS.experimentID{I}; N = length(experimentID);
        experimentID_D = D.experimentID{I};
        experimentID{N+1} = experimentID_D{1};
        DS.experimentID{I} = experimentID;
    end
end
% 
D = DS;