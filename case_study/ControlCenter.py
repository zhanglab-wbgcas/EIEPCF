import pandas as pd
import numpy as np
from EIEPCF import EIEPCF
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

data_name = "Arab_coldResistance"
print('--------------------' + data_name + '--------------------')

# ----------------------------------------------Data preprocessing------------------------------------------------

expression_data_raw = pd.read_csv('./GSE162507_all_FPKM.csv', sep=',', index_col=0)
gene_name_all = expression_data_raw.index.values

gold_standard = pd.read_csv('./GoldStandard.csv')
gold_standard = gold_standard.values

regulator_name_raw = gold_standard[:, 0]
gene_name_raw = gold_standard[:, 1]
regulator_name_raw = np.unique(regulator_name_raw)
gene_name_raw = np.unique(gene_name_raw)

regulator_name = []
for i in regulator_name_raw:
    if (gene_name_all.__contains__(i)):
        regulator_name.append(i)

gene_name = []
for i in gene_name_raw:
    if (gene_name_all.__contains__(i)):
        gene_name.append(i)

all_names = []
all_names.extend(regulator_name)
all_names.extend(gene_name)
all_names = np.unique(all_names)

print("The number of regulators:", len(regulator_name))
print("The number of genes:", len(all_names))

expression_data = expression_data_raw.loc[all_names, :]
expression_data.to_csv('./log/' + data_name + '_ExpressionData.csv')
expression_data = expression_data.values.T
expression_data = np.log1p(expression_data)

gold_matrix = np.zeros((len(regulator_name), len(all_names)), dtype=int)
gold_matrix = pd.DataFrame(gold_matrix)
gold_matrix.columns = all_names
gold_matrix.index = regulator_name

for i in range(len(gold_standard)):
    temp = gold_standard[i]
    if (gene_name.__contains__(temp[0]) and gene_name.__contains__(temp[1])):
        gold_matrix[temp[1]][temp[0]] = 1

gold_matrix.to_csv('./log/' + data_name + '_GoldStandard.csv')

# ----------------------------------------------Calculate and save the regulation matrix-------------------------------

regu_matrix = EIEPCF(expression_data, all_names.tolist(), regulator_name)

regu_matrix = pd.DataFrame(regu_matrix)
regu_matrix.columns = all_names
regu_matrix.index = regulator_name
regu_matrix.to_csv('./log/' + data_name + '_EIEPCF_Matrix.csv')

# ---------------------------------------------Draw venn diagram-------------------------------------------------------

predicted = []
gold_value = []

for i in all_names.tolist():
    for j in regulator_name:
        predicted.append(regu_matrix[i][j])
        gold_value.append(gold_matrix[i][j])

predicted = np.asarray(predicted)
gold_value = np.asarray(gold_value)

abs_predicted = abs(predicted)

abs_predicted2 = []
c = 0
que = 0.8
for i in abs_predicted:
    if (i > que):
        abs_predicted2.append(1)
        c = c + 1
    else:
        abs_predicted2.append(0)

regu_matrix_bi = pd.DataFrame(regu_matrix)
regu_matrix_bi.columns = all_names
regu_matrix_bi.index = regulator_name
all_nets = []
all_nets1 = []
id2name = pd.read_csv("id_geneName.csv", index_col=0)

regu_matrix_bi_all = []
for i in all_names.tolist():
    for j in regulator_name:
        l = 0
        if (abs(regu_matrix[i][j]) > que):
            s = j + "," + i
            regu_matrix_bi_all.append(s)
            regu_matrix_bi[i][j] = 1
            l = 1
            jn = id2name.loc[j].values[0]
            inn = id2name.loc[i].values[0]
            net1 = (j, i, "tf", "target")
            all_nets1.append(net1)
            net = (jn, inn, "tf", "target")
            all_nets.append(net)
        else:
            regu_matrix_bi[i][j] = 0
            l = 0

pd.DataFrame(all_nets1).to_csv('./log/' + data_name + "_all_nets_end.csv")
pd.DataFrame(all_nets).to_csv('./log/' + data_name + "_all_nets_name_end.csv")
regu_matrix_bi.to_csv('./log/' + data_name + '_EIEPCF_Matrix_bi.csv')

gold_standard_all = []
for i in gold_standard:
    s = i[0] + "," + i[1]
    gold_standard_all.append(s)

set1 = set(gold_standard_all)
set2 = set(regu_matrix_bi_all)
intersection = set1 & set2
print("The number of intersections:", len(intersection))
venn2(subsets=[set1, set2], set_labels=('Reported ', 'Prediction'))
plt.show()

# ---------------------------------------------------------------------------------------------------------------------
