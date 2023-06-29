# EIEPCF
1. Construct TF-TF gene regulatory network
TF-TF gene regulatory network means that all genes constructing regulatory network can regulate other genes, that is, all genes have the function of TF. Example code for constructing TF-TF gene regulation network is as follows:
import pandas as pd
from EIEPCF import EIEPCF

# Read gene expression data
expression_data = pd.read_csv('./Data/GeneExpressionData.csv', sep=',', index_col=0)
# Get gene names and convert gene names to a list
gene_name = expression_data.index.values.tolist()
# Assign a value to regulator names
regulator_name = gene_name
# Get the gene expression matrix (rows are samples, columns are genes)
expression_data = expression_data.values.T

# The regulatory relationships between genes were inferred by EIEPCF method
# The first way
regu_matrix = EIEPCF(expression_data)
# The second way
regu_matrix = EIEPCF(expression_data, gene_name=gene_name)
# The third way
regu_matrix = EIEPCF(expression_data, gene_name=gene_name, regulator=regulator_name)

# Save the regulation matrix
regu_matrix = pd.DataFrame(regu_matrix)
regu_matrix.columns = gene_name
regu_matrix.index = regulator_name
regu_matrix.to_csv('./log/TF_TF_ReguMatrix_EIEPCF.csv')

2. Construct TF-Gene gene regulatory network
TF-Gene gene regulatory network means that some genes constructing regulatory network can regulate other genes, that is, have the function of TF; Some can only be regulated by other genes, but cannot regulate other genes, that is, they do not have the function of TF. Example code for constructing TF-Gene gene regulation network is as follows:
import pandas as pd
from EIEPCF import EIEPCF

# Read gene expression data
expression_data = pd.read_csv('./Data/GeneExpressionData.csv', sep=',', index_col=0)
# Get gene names and convert gene names to a list
gene_name = expression_data.index.values.tolist()
# Read regulator names (all regulator names must be included in gene names)
regulator_name = pd.read_csv(' ./Data/RegulatorName.csv', sep=',', header=None)
# Convert regulator names to a list
regulator_name = regulator_name.values
regulator_name = regulator_name[0].tolist()
# Get the gene expression matrix (rows are samples, columns are genes)
expression_data = expression_data.values.T

# The regulatory relationships between genes were inferred by EIEPCF method
regu_matrix = EIEPCF(expression_data, gene_name=gene_name, regulator=regulator_name)

# Save the regulation matrix
regu_matrix = pd.DataFrame(regu_matrix)
regu_matrix.columns = gene_name
regu_matrix.index = regulator_name
regu_matrix.to_csv('./log/TF_Gene_ReguMatrix_EIEPCF.csv')
