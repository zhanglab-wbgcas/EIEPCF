import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
from GENIE3 import GENIE3

data_name = "simData_size100"
print('--------------------' + data_name + '--------------------')

# ----------------------------------------------Read data--------------------------------------------------------------

# Read the gene expression data
gene_expression = pd.read_csv('E:/Code/MATLAB/Simulate_Data/log/simData/' + data_name + '_Y.csv', sep=',', header=None)
gene_expression = gene_expression.values.T
# Read the regulator expression data
regulator_expression = pd.read_csv('E:/Code/MATLAB/Simulate_Data/log/simData/' + data_name + '_X.csv', sep=',', header=None)
regulator_expression = regulator_expression.values.T
# Read the gold standard, the rows are the genes, the columns are the regulators
gold_matrix = pd.read_csv('E:/Code/MATLAB/Simulate_Data/log/simData/' + data_name + '_Gold.csv', sep=',', header=None)
gold_matrix = gold_matrix.values

# Gets the number of genes
gene_number = gene_expression.shape[1]
# Gets the number of regulators
regulator_number = regulator_expression.shape[1]

# ----------------------------------------------Gets the regulatory matrix--------------------------------------------

expression_data = np.hstack((gene_expression, regulator_expression))
gene_name = list(range(gene_number))
regulator_name = list(range(gene_number, expression_data.shape[1]))
all_name = list(range(expression_data.shape[1]))

regu_matrix = GENIE3(expression_data, gene_names=all_name, regulators=regulator_name)
regu_matrix = regu_matrix[gene_number:, 0:gene_number]

regu_matrix = regu_matrix.T

# ----------------------------------------------Save regulatory matrix------------------------------------------------

regu_matrix_save = pd.DataFrame(regu_matrix)
regu_matrix_save.to_csv('./log/ReguMatrix/' + data_name + '_RF_Matrix.csv')

# ----------------------------------------------Perform performance evaluation-----------------------------------------

predicted = []
gold_value = []
for i in range(gene_number):
    for j in range(regulator_number):
        predicted.append(regu_matrix[i, j])
        gold_value.append(gold_matrix[i, j])

predicted = np.asarray(predicted)
gold_value = np.asarray(gold_value)

fpr, tpr, threshold = roc_curve(gold_value, abs(predicted), pos_label=1)
AUROC = auc(fpr, tpr)
print("AUROC:", AUROC)

# ------------------------------------------------------------------------------------------------------------
