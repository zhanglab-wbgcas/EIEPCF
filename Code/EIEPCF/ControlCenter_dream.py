import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
from EIEPCF import EIEPCF

data_name = "Dream3_E2_50"
print('--------------------' + data_name + '--------------------')

# ----------------------------------------------Gets the gene regulation matrix---------------------------------------

expression_data = pd.read_csv('E:/Data/Dream3Data/InSilicoSize50-Ecoli2-null-mutants.tsv', sep='\t', index_col=0)
expression_data = expression_data.drop(['wt'])

gene_name = expression_data.columns.values.tolist()
regulator_name = gene_name
expression_data = expression_data.values

regu_matrix = EIEPCF(expression_data)

# ----------------------------------------------Save the regulation matrix---------------------------------------------

regu_matrix = pd.DataFrame(regu_matrix)
regu_matrix.columns = gene_name
regu_matrix.index = regulator_name
regu_matrix.to_csv('./log/EIEPCF_ReguMatrix/' + data_name + '_EIEPCF_Matrix.csv')

# ---------------------------------------------Perform performance evaluation------------------------------------------

gold_standard = pd.read_csv('E:/Data/Dream3Data/DREAM3GoldStandard_InSilicoSize50_Ecoli2.txt', sep='	', header=None)
gold_standard = gold_standard.values

gold_value = gold_standard[:, 2]
gold_value = np.asarray(gold_value, dtype=int)

predicted = []
# 获取预测结果
for i in range(len(gold_value)):
    predicted.append(regu_matrix[gold_standard[i, 1]][gold_standard[i, 0]])
predicted = np.asarray(predicted)

# ROC曲线下面积（AUC）
fpr, tpr, threshold = roc_curve(gold_value, abs(predicted), pos_label=1)
AUROC = auc(fpr, tpr)
print("AUROC:", AUROC)

