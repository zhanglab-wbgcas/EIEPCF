function [J0,X,Y]=simdata(n_gene,m_sample,p_TF,noise,rate)
    % Method:
    % 模拟表达数据和基准网络
    % Input:
    % n_gene:基因的尺寸
    % m_sample:样本的尺寸
    % p_TF:调控因子的个数
    % noise:添加噪声时，噪声的方差
    % rate:生成网络时，网络边在全连接网络中所占的比例
   
    % Output:
    % J0:生成的基准网络，行为靶基因，列为调控因子
    % X:模拟的调控因子的表达量，行为基因，列为样本
    % Y:模拟的靶基因的表达量，行为基因，列为样本
    
    % Author: PengHuixiang
    % Version data: 2023-02-15
    
    %%
    % rand()是0-1的均匀随机分布，randn()是均值为0方差为1的正态分布
    
    % 为随机生成器设置随机种子
    rng(4)
    % 生成一个随机分布的矩阵
    J0=rand(n_gene,p_TF);
    % 对矩阵中所有的元素进行排序
    A = sort(J0(:)); 
    % 根据设置的比率获取调控矩阵中的0、1值，即挑选具有调控关系的边
    t = A(size(J0,1)*size(J0,2)*rate);
    J0 = (J0 < t);
    
    % 生成均值为1方差为1的高斯分布值X，即调控因子的表达量
    X=1+sqrt(1)*randn(p_TF,m_sample);
    
    % 生成无噪声的靶基因的表达量
    Y=J0*X;
    % 生成均值为0方差为noise的高斯噪声
    a=0+sqrt(noise)*randn(n_gene,m_sample);
    % 生成由噪声的靶基因的表达量
    Y=Y+a;
    
end





