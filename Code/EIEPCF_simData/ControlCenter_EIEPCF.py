import pandas as pd
import numpy as np
from sklearn.linear_model import Lasso
from sklearn.metrics import roc_curve, auc

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

# ----------------------------------------------Gets the regulatory matrix---------------------------------------------

regu_matrix = []

for i in range(gene_number):
    temp = []
    for j in range(regulator_number):
        Y = gene_expression[:, i]
        T = regulator_expression[:, j]
        X = np.delete(regulator_expression, [j], axis=1)

        model_Y = Lasso()
        model_T = Lasso()

        model_Y.fit(X, Y)
        model_T.fit(X, T)

        Y1 = model_Y.predict(X)
        T1 = model_T.predict(X)

        Y2 = Y - Y1
        T2 = T - T1

        ate_ = np.correlate(Y2, T2)
        ate_ = ate_[0]

        print(str(j) + "_" + str(i) + ":", ate_)

        temp.append(ate_)

    regu_matrix.append(temp)

regu_matrix = np.asarray(regu_matrix)

# ----------------------------------------------Save regulatory matrix------------------------------------------------

regu_matrix_save = pd.DataFrame(regu_matrix)
regu_matrix_save.to_csv('./log/ReguMatrix/' + data_name + '_EIEPCF_Matrix.csv')

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
