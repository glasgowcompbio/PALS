import copy
import random


class Method(object):
    def __init__(self, data_source, seed=None, preprocessors=None):
        if seed is None:
            random.seed()
        else:
            random.seed(seed)

        self.data_source = copy.deepcopy(data_source)
        if preprocessors is None:
            self.preprocessors = self._create_default_preprocessors()
        else:
            self.preprocessors = preprocessors

    def get_results(self, preprocess=True):
        raise NotImplementedError()

    def get_pathway_df(self, standardize=True): # to support old method calls
        return self.get_results(preprocess=standardize)

    def _create_default_preprocessors(self):
        return []

    def _get_measurement_df(self, preprocess):
        df = self.data_source.get_measurements()
        if preprocess:
            for p in self.preprocessors:
                df = p.process(df)
        return df