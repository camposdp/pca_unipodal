%% Cálculo de PCA de ângulos articulares
% Daniel Prado de Campos - 02/2022 
% UTFPR-PR/PUC-PR 

%% Limpar workspace
clear
clc
close all

%% Definições iniciais


angle_set = [1 5 ... %hip (dois lados) 
             2 6 ... %knee (dois lados)
             4 7 ... %ankle (dois lados)
             10 ...  % pelvis(direito)
             13 ...  %neck(direito)
             15 ...  %spine (direito)
             23]; %thorax (direito)
             %25]; % head (direito) 

fs= 200; %frequência de amostragem



%Cortar dados a partir de qual %?
%e.g. de 10% a 90%: xi = 0.1, xf = 0.9
% xi = 0.1;  
% xf = 0.9; 

%em segundos
ti = 15;  
tf = 50; 


%% Abrir dados

path = 'D:\PESQUISA\UCM\MATLAB\_importarPyCGM_compararVicon\pilotoMarcelo\Marcelo\19_11_21\';
file = 'UNIP_DIR_1';
ext = '.c3d';

[AQ, byteOrder, storageFormat] = btkReadAcquisition([path file ext]);
Markers=btkGetMarkers(AQ);
Angles=btkGetAngles(AQ);

%% Cálculo da matriz normalizada

angle_names_all=fieldnames(Angles);
for i = angle_set
selected_angles.(angle_names_all{i}) = Angles.(angle_names_all{i}); 
end
angle_names=fieldnames(selected_angles);
angle_matrix_raw = cell2mat(struct2cell(selected_angles)');
angle_matrix_raw( ~any(angle_matrix_raw,2), : ) = []; 
N = length(angle_matrix_raw);
M = length(angle_names);
%angle_matrix=angle_matrix(floor(N*xi):floor(N*xf),:);


%% user input
fprintf('\n')
fprintf('0 - Todas')
for k = 1:M
    fprintf('\n')
    fprintf(num2str(k))
    fprintf(' - ')
    fprintf(angle_names{k})
end
fprintf('\n')
fprintf('\n')
prompt1 = 'Analisar qual ângulo? Digite o número: ';
x1 = input(prompt1);
fprintf('\n')
prompt2 = 'Qual componente analisar? Digite o número: ';
x2 = input(prompt2);
fprintf('\n')

index = x1;
PC = x2;


%% Projeto do filtro
N=2; %ordem do filtro
fNyquist = fs/2;

fCorte = 10;
fCorte_normalizada = fCorte/fNyquist; 

[B,A] = butter(N,fCorte_normalizada,'low'); % IIR

s = size(angle_matrix_raw);
for i=1:s(2)
    angle_matrix(:,i) = filtfilt(B,A,angle_matrix_raw(:,i));
end

%%


angle_matrix=angle_matrix(floor(ti*fs):floor(tf*fs),:);
angle_mean = mean(angle_matrix);
d_norm=vecnorm(angle_matrix,2,2); % é escalar ou vetor?
%Se escalar:
%d_norm=norm(mean(angle_matrix));

angle_norm = (angle_matrix-angle_mean)./d_norm;
%angle_norm = (bsxfun(@minus, angle_matrix, angle_mean))./d_norm;

%% Cálculo da PCA e da reconstrução


 %Qual componente analisar?


[coeff,score,latent,tsquared,explained,mu] = pca(angle_norm);
PM = angle_mean + d_norm.*score(:,PC)*coeff(:,PC)';
%PM = angle_mean + d_norm.*score(:,1:2)*coeff(:,1:2)';

%% Gráfico comparando original vs projeção
t = 0:1/fs:(length(angle_matrix)-1)/fs;
if x1 > 0
ix = (index-1)*3+1;


figure
subplot(3,1,1)
hold all
%grid on
plot(t,angle_matrix(:,ix),'-','color',[0.8 0.8 0.8],'LineWidth',1.2)
plot(t,PM(:,ix),'r','LineWidth',1.2)
xlabel('Tempo (s)')
ylabel('Ângulo')
s = ['PC = ' num2str(PC) ' - ' angle_names{index} ' - X '];
legend('Original','Componente','Location','Best')
title(s);

subplot(3,1,2)
hold all
%grid on
plot(t,angle_matrix(:,ix+1),'-','color',[0.8 0.8 0.8],'LineWidth',1.2)
plot(t,PM(:,ix+1),'g','LineWidth',1.2)
xlabel('Tempo (s)')
ylabel('Ângulo')
s = ['PC = ' num2str(PC) ' - ' angle_names{index} ' - Y '];
legend('Original','Componente','Location','Best')
title(s);

subplot(3,1,3)
hold all
%grid on
plot(t,angle_matrix(:,ix+2),'-','color',[0.8 0.8 0.8],'LineWidth',1.2)
plot(t,PM(:,ix+2),'b','LineWidth',1.2)
xlabel('Tempo (s)')
ylabel('Ângulo')
s = ['PC = ' num2str(PC) ' - ' angle_names{index} ' - Z '];
legend('Original','Componente','Location','Best')
title(s);

%%
else
f2 = figure
tile = tiledlayout(3,M,'TileSpacing','tight');
xlabel(tile,'Tempo (s)')
ylabel(tile,'Ângulo')
title_string = ['PC = ' num2str(PC) ', Explicado: ' num2str(explained(PC),'%.2f') '%'];
title(tile,title_string)
fSize = 6;
for j = 1:M
index = j;   
ix = (j-1)*3+1;   
    
%subplot(M,3,ix)
nexttile(j)
hold all
%grid on
plot(t,angle_matrix(:,ix),'-','color',[0.8 0.8 0.8],'LineWidth',1.2)
plot(t,PM(:,ix),'r','LineWidth',1.2)

%xlabel('Tempo (s)')
%ylabel('Ângulo')
%s = ['PC = ' num2str(PC) ' - ' angle_names{index} ' - X '];
%legend('Original','Componente')
s = [angle_names{index}];
title(s);
ax = gca;
ax.FontSize = fSize; 

%subplot(M,3,ix+1)
nexttile(M+j)
hold all
%grid on
plot(t,angle_matrix(:,ix+1),'-','color',[0.8 0.8 0.8],'LineWidth',1.2)
plot(t,PM(:,ix+1),'g','LineWidth',1.2)
%xlabel('Tempo (s)')
%ylabel('Ângulo')
s = ['PC = ' num2str(PC) ' - ' angle_names{index} ' - Y '];
%legend('Original','Componente')
%title(s);
ax = gca;
ax.FontSize = fSize; 

%subplot(M,3,ix+2)
nexttile(2*M+j)
hold all
%grid on
plot(t,angle_matrix(:,ix+2),'-','color',[0.8 0.8 0.8],'LineWidth',1.2)
plot(t,PM(:,ix+2),'b','LineWidth',1.2)
%xlabel('Tempo (s)')
%ylabel('Ângulo')
s = ['PC = ' num2str(PC) ' - ' angle_names{index} ' - Z '];
%legend('Original','Componente')
%title(s);
ax = gca;
ax.FontSize = fSize;    

f2.Position = [100 80 1200 600];

[c1,p1] = corrcoef(angle_matrix(:,ix),PM(:,ix));
[c2,p2] = corrcoef(angle_matrix(:,ix+1),PM(:,ix+1));
[c3,p3] = corrcoef(angle_matrix(:,ix+2),PM(:,ix+2));
R(1,j)= c1(2,1);
R(2,j)= c2(2,1);
R(3,j)= c3(2,1);

P(1,j)= p1(2,1);
P(2,j)= p2(2,1);
P(3,j)= p3(2,1);


end
set(gcf, 'Renderer', 'Painters');
print([file '_GRAFICOS_C' num2str(PC)],'-depsc2');
CROSS_CORR = R;

%% HEATMAP

%Mapas de cor
%https://www.mathworks.com/help/matlab/ref/colormap.html
% Para inverter a escala usar: flipud(hot)

f1=figure
eixo_y = {'x','y','z'};
h=heatmap(angle_names,eixo_y,R,'Colormap',turbo)
ylabel('Eixo')
xlabel('Ângulo')

f1.Position = [100 100 1000 350];
h.CellLabelFormat = '%.2f';
print([file '_HEATMAP_C' num2str(PC)],'-depsc2');

end

