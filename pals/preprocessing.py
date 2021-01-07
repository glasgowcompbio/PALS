import warnings

import numpy as np
import pandas as pd
from loguru import logger
from sklearn import preprocessing
from sklearn.impute import KNNImputer


class Preprocessing(object):
    def process(self):
        raise NotImplementedError()


class ZeroAndNegativeReplace(Preprocessing):
    def process(self, df):
        logger.debug('Replacing negative and zero values with NaN')
        df[df < 0] = 0  # prevent negative values in df
        df[df == 0.0] = None  # replace 0.0 with NaN for easier operations ahead
        return df


class MinValueImputation(Preprocessing):
    def __init__(self, groups, min_replace):
        self.groups = groups
        self.min_replace = min_replace

    def process(self, df):
        logger.debug('Performing min-value imputation')
        # if all intensities in a (factor) group are zero, a min value is set.
        for group_name, samples in self.groups.items():
            # If all zero in group then replace with minimum
            df.loc[df.loc[:, samples].isnull().all(axis=1), samples] = self.min_replace
        return df


class RowAverageImputation(Preprocessing):
    def __init__(self, groups):
        self.groups = groups

    def process(self, df):
        logger.debug('Performing row average imputation')
        for group_name, samples in self.groups.items():
            # if group values not all zeros, replace the zeros with mean of group
            group_df = df.loc[:, samples]
            df.loc[:, samples] = group_df.mask(group_df.isnull(), group_df.mean(axis=1), axis=0)
        return df


class KNNImputation(Preprocessing):
    def __init__(self, groups, K=5):
        self.groups = groups
        self.K = K

    def process(self, df):
        logger.debug('Performing K-Nearest Neighbour imputation')
        for group_name, samples in self.groups.items():
            group_df = df.loc[:, samples]
            # if group values not all zeros, perform KNN imputation
            if np.sum(group_df.isnull().values) > 0:
                imputer = KNNImputer(n_neighbors=self.K)
                imputed_df = pd.DataFrame(imputer.fit_transform(group_df), index=group_df.index,
                                          columns=group_df.columns)
                df.loc[:, samples] = imputed_df
        return df


class LogNormalisation(Preprocessing):
    def process(self, df):
        logger.debug('Applying log normalisation')
        return np.log(df)


class ZScoreNormalisation(Preprocessing):
    def process(self, df):
        # standardize the data across the samples (zero mean and unit variance))
        logger.debug('Scaling the data across the sample: zero mean and unit variance')

        # https://stackoverflow.com/questions/15933741/how-do-i-catch-a-numpy-warning-like-its-an-exception-not-just-for-testing
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                # for some reason, calling preprocessing.scale() below sometimes throws the following warning:
                # > "UserWarning: Numerical issues were encountered when scaling the data and might not be solved.
                # > The standard deviation of the data is probably very close to 0."
                #
                # This seems to be fixed if we use the StandardScaler instead, which does the same thing
                # see https://stackoverflow.com/questions/51741605/standardize-dataset-containing-too-large-values

                # scaled_data = np.log(np.array(df))
                # scaled_data = preprocessing.scale(scaled_data, axis=1)
                # sample_names = df.columns
                # df[sample_names] = scaled_data

                scaler = preprocessing.StandardScaler()
                scaled_df = scaler.fit_transform(df.transpose())  # transpose into the right shape for StandardScaler
                df = pd.DataFrame(scaled_df.transpose(), columns=df.columns,
                                  index=df.index)  # return the original shape
                return df
            except UserWarning as e:
                raise (e)
