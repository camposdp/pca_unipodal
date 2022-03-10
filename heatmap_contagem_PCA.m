%%% Análise dos dados da PCA processados
% Daniel Prado de Campos
% UTFPR/PUC-PR
% 10/03/2022


%% Inicializar

clc
close all
clear


load dados_pca PCA*

lado = 'E'; %D/E
dim = 'Z'; %X/Y/Z
th = 0.8; %Limiar do cross-correlation
%% Threshold 


for comp = 1:8

str_var = ['PCA' num2str(comp) '_' lado '_' dim];

X = eval(str_var);
count_th(comp,:) = sum(X>th);
end

eixo_y = {'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'};

angle_names = {'Hip_L','Hip_R','Knee_L','Knee_R','Ankle_L','Ankle_R','Pelvis','Neck','Spine','Thorax'};
figure
heatmap(angle_names,eixo_y,count_th);

%% Greater

for comp = 1:8

str_var = ['PCA' num2str(comp) '_' lado '_' dim];

X = eval(str_var);
[~,idx]=max(X');
[a,b]=hist(idx,unique(idx));
for i = 1:length(a)
count_max(comp,b) = a; 
end

end

eixo_y = {'PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'};
angle_names = {'Hip_L','Hip_R','Knee_L','Knee_R','Ankle_L','Ankle_R','Pelvis','Neck','Spine','Thorax'};
figure
heatmap(angle_names,eixo_y,count_max);

