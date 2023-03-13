% Use: ControlCenter
% Version data: 2023-02-15
% Author: PengHuixiang

% Empty the variable
clear
% Clear the command line window
clc
% Modify the path as "E:/Code/MATLAB/Simulate_Data"
cd E:/Code/MATLAB/Simulate_Data
% Show the path
disp(pwd)
% Add the directories "Simulate_Data" into MATLABPATH
path(path,'E:/Code/MATLAB/Simulate_Data')


%% Size10

% % The network size is 10
% data_name = 'simData_size10';
% % Set parameters
% genes_number = 10;
% sample_number = 5;  
% TFs_number = 10; 
% noise = 0.1;    
% rate  = 0.1; 

% % The network size is 50
% data_name = 'simData_size50';
% % Set parameters
% genes_number = 10;
% sample_number = 8;  
% TFs_number = 50; 
% load 50sim1261649_10.mat
% Gold = double(J0);

% % The network size is 100
data_name = 'simData_size100';
% Set parameters
genes_number = 10;
sample_number = 10;  
TFs_number = 100; 
noise = 0.1;    
rate  = 0.05; 

% % The network size is 1000
% data_name = 'simData_size1000';
% % Set parameters
% genes_number = 10;
% sample_number = 20;  
% TFs_number = 1000; 
% noise = 0.1;    
% rate  = 0.002;


% Simulated expression data and benchmark networks
[Gold,X,Y]=simdata(genes_number,sample_number,TFs_number,noise,rate);
Gold = double(Gold);

% Save
filename = strcat('./log/simData/',data_name,'_Gold.csv');
csvwrite(filename, Gold)
filename = strcat('./log/simData/',data_name,'_X.csv');
csvwrite(filename, X)
filename = strcat('./log/simData/',data_name,'_Y.csv');
csvwrite(filename, Y)


%% LP

method_name = 'LP';

% Obtain a regulatory matrix with rows for genes and columns for transcription factors
lamda = 1;
Regu_Matrix = zeros(genes_number,TFs_number);
for i=1:genes_number
    y = Y(i,:);
    Regu_Matrix(i,:) = LP_TGN(y,X,lamda);
end

filename = strcat('./log/ReguMatrix/',data_name,'_',method_name,'_Matrix.csv');
csvwrite(filename, Regu_Matrix)

Regu_Matrix = abs(Regu_Matrix);

gold_edge = zeros(genes_number*TFs_number,1);
pre_value = zeros(genes_number*TFs_number,1);
flag = 0;
for i = 1:genes_number
    for j = 1:TFs_number
        flag = flag + 1;
        gold_edge(flag,1) = Gold(i,j);
        pre_value(flag,1) = Regu_Matrix(i,j);
    end
end

roc = ROC([gold_edge, pre_value], 0);
fprintf('AUROC: %.4f \n', roc.auc);


%% LASSO

method_name = 'LASSO';

% Obtain a regulatory matrix with rows for genes and columns for transcription factors
lamda = 0.5;
Regu_Matrix = zeros(genes_number,TFs_number);
for i=1:genes_number
    disp(i);
    y = Y(i,:);
    temp=larsen(X', y', lamda);
    Regu_Matrix(i,:) = temp(size(temp,1),:);
end

filename = strcat('./log/ReguMatrix/',data_name,'_',method_name,'_Matrix.csv');
csvwrite(filename, Regu_Matrix)

Regu_Matrix = abs(Regu_Matrix);

gold_edge = zeros(genes_number*TFs_number,1);
pre_value = zeros(genes_number*TFs_number,1);
flag = 0;
for i = 1:genes_number
    for j = 1:TFs_number
        flag = flag + 1;
        gold_edge(flag,1) = Gold(i,j);
        pre_value(flag,1) = Regu_Matrix(i,j);
    end
end

roc = ROC([gold_edge, pre_value], 0);
fprintf('AUROC: %.4f \n', roc.auc);


%% PCA-CMI

method_name = 'PCA-CMI';

order = 0;
Regu_Matrix = zeros(genes_number,TFs_number);

for i = 1:genes_number
    for j = 1:TFs_number
        Regu_Matrix(i,j)=cmi(Y(i,:),X(j,:));
    end
end

filename = strcat('./log/ReguMatrix/',data_name,'_',method_name,'_Matrix.csv');
csvwrite(filename, Regu_Matrix)

gold_edge = zeros(genes_number*TFs_number,1);
pre_value = zeros(genes_number*TFs_number,1);
flag = 0;
for i = 1:genes_number
    for j = 1:TFs_number
        flag = flag + 1;
        gold_edge(flag,1) = Gold(i,j);
        pre_value(flag,1) = Regu_Matrix(i,j);
    end
end

roc = ROC([gold_edge, pre_value], 0);
fprintf('AUROC: %.4f \n', roc.auc);


%% NARROMI

method_name = 'NARROMI';

lamda =  1; 
alpha = 0.05;
beta = 0.05;
t = 0.6; 

Regu_Matrix = zeros(genes_number,TFs_number);
for i=1:genes_number
    y = Y(i,:);
    [~,temp,~]=narromi(y', X', lamda, alpha, beta, t);
    Regu_Matrix(i,:) = temp;
end

filename = strcat('./log/ReguMatrix/',data_name,'_',method_name,'_Matrix.csv');
csvwrite(filename, Regu_Matrix)

gold_edge = zeros(genes_number*TFs_number,1);
pre_value = zeros(genes_number*TFs_number,1);
flag = 0;
for i = 1:genes_number
    for j = 1:TFs_number
        flag = flag + 1;
        gold_edge(flag,1) = Gold(i,j);
        pre_value(flag,1) = Regu_Matrix(i,j);
    end
end

roc = ROC([gold_edge, pre_value], 0);
fprintf('AUROC: %.4f \n', roc.auc);


%%











