import numpy as np
from numpy import ndarray
from sklearn.linear_model import Lasso
import time


def EIEPCF(expr_data, gene_name=None, regulator='all'):
    '''

    Calculate correlations based on causal inference

    :param expr_data:numpy array
        Gene expression data, row represents sample, column represents gene

    :param gene_name:list of strings, optional
        Name of gene
        default: None

    :param regulator:list of strings, optional
        Name of regulatory factor
        default: 'all'

    :param thread_number:positive integer, optional
        Number of threads used for parallel computing
        default: 1

    :return:


    '''

    time_start = time.time()

    # Check input arguments
    if not isinstance(expr_data, ndarray):
        raise ValueError(
            'expr_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')

    gene_number = expr_data.shape[1]

    if gene_name is not None:
        if not isinstance(gene_name, (list, tuple)):
            raise ValueError('input argument gene_name must be a list of gene names')
        elif len(gene_name) != gene_number:
            raise ValueError(
                'input argument gene_name must be a list of length p, where p is the number of columns/genes in the expr_data')

    if regulator != 'all':
        if not isinstance(regulator, (list, tuple)):
            raise ValueError('input argument regulator must be a list of gene names')

        if gene_name is None:
            raise ValueError('the gene name must be specified (in input argument gene_name)')
        else:
            sIntersection = set(gene_name).intersection(set(regulator))
            if not sIntersection:
                raise ValueError('the genes must contain at least one candidate regulator')

    # Get the indices of the candidate regulators
    if regulator == 'all':
        regulator_idx = list(range(gene_number))
    else:
        regulator_idx = [i for i, gene in enumerate(gene_name) if gene in regulator]

    regulator_number = len(regulator_idx)

    regu_matrix = []

    for i in range(gene_number):
        temp = []
        for j in range(regulator_number):
            Y = expr_data[:, i]
            T = expr_data[:, regulator_idx[j]]

            k = regulator_idx[j]
            X = np.delete(expr_data, [i, k], axis=1)

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

            if i == k:
                ate_ = 0
            temp.append(ate_)

        regu_matrix.append(temp)

    regu_matrix = np.asarray(regu_matrix).T

    time_end = time.time()
    print("Elapsed time: %.2f seconds" % (time_end - time_start))

    return regu_matrix
