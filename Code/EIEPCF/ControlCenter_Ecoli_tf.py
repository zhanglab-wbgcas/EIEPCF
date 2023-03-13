import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
from EIEPCF import EIEPCF

data_name = "Ecoli_tf"
print('--------------------' + data_name + '--------------------')

# ----------------------------------------------Data preprocessing------------------------------------------------

expression_data_raw = pd.read_csv('E:/Data/Ecoli_Real/avg_E_coli_v4_Build_6_exps466probes7459.tab', sep='\t')

gene_name_raw = list(expression_data_raw['E_coli_v4_Build_6:allProbes'])

gene_name = []
for i in gene_name_raw:
    temp = i.split('_')[0]
    temp = temp.lower()
    gene_name.append(temp)

gold_standard_raw = pd.read_csv('E:/Data/Ecoli_Real/network_tf_tf.txt', sep='\t', header=None)
gold_standard_raw = gold_standard_raw.values

gold_standard = []
for i in range(len(gold_standard_raw)):
    temp = gold_standard_raw[i]
    temp[0] = temp[0].lower()
    temp[1] = temp[1].lower()
    if temp[0] != temp[1]:
        if temp[0] in gene_name:
            if temp[1] in gene_name:
                gold_standard.append(temp)
gold_standard = np.asarray(gold_standard)

regulator_name = []
regulator_name.extend(gold_standard[:, 0])
regulator_name.extend(gold_standard[:, 1])
regulator_name = np.unique(regulator_name)
regulator_name = regulator_name.tolist()

expression_data_raw.index = gene_name
expression_data_raw = expression_data_raw.drop(columns=['E_coli_v4_Build_6:allProbes'], axis=1)
expression_data = expression_data_raw.loc[regulator_name, :]
expression_data.to_csv('./log/Ecoli_Real/' + data_name + '_ExpressionData.csv')

expression_data = expression_data.values.T

# ----------------------------------------------Calculate and save the regulation matrix-------------------------------

regu_matrix = EIEPCF(expression_data)

regu_matrix = pd.DataFrame(regu_matrix)
regu_matrix.columns = regulator_name
regu_matrix.index = regulator_name
regu_matrix.to_csv('./log/Ecoli_Real/' + data_name + '_EIEPCF_Matrix.csv')

# ---------------------------------------------Perform performance evaluation----------------------------------------

gold_matrix = np.zeros((len(regulator_name), len(regulator_name)), dtype=int)
gold_matrix = pd.DataFrame(gold_matrix)
gold_matrix.columns = regulator_name
gold_matrix.index = regulator_name
for i in range(len(gold_standard)):
    temp = gold_standard[i]
    gold_matrix[temp[1]][temp[0]] = 1
gold_matrix.to_csv('./log/Ecoli_Real/' + data_name + '_GoldStandard.csv')

predicted = []
gold_value = []
for i in range(len(regulator_name)):
    for j in range(len(regulator_name)):
        predicted.append(regu_matrix[regulator_name[j]][regulator_name[i]])
        gold_value.append(gold_matrix[regulator_name[j]][regulator_name[i]])

predicted = np.asarray(predicted)
gold_value = np.asarray(gold_value)

fpr, tpr, threshold = roc_curve(gold_value, abs(predicted), pos_label=1)
AUC = auc(fpr, tpr)
print("AUC:", AUC)

# ---------------------------------------------------------------------------------------
