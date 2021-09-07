from mummichog import get_user_data
from mummichog.functional_analysis import InputUserData, DataMeetModel, PathwayAnalysis
from mummichog.models import metabolicNetwork, metabolicModels
from mummichog.reporting import LocalExporting
import pandas as pd
from pals.base import Method
import os
import numpy as np
import matplotlib.pyplot as plt

class MummichogDataExporting(LocalExporting):

    def __init__(self, pathway_analysis):
        self.pathway_analysis = pathway_analysis
        self.pathway_enrich_df = pd.DataFrame(columns=['pathway', 'overlap_size', 'pathway_size', 'p-value',
                                                       'overlap_Empirical', 'overlap_features (id)',
                                                       'overlap_features (name)'])

    def construct_pathway_enrich_df(self):

        resultstr = [['pathway', 'overlap_size', 'pathway_size', 'p-value',
                                                       'overlap_Empirical', 'overlap_features (id)',
                                                       'overlap_features (name)']]

        for P in self.pathway_analysis.resultListOfPathways:

            comp = P.overlap_EmpiricalCompounds
            EIDs = []
            cpds = []

            for EID_set in P.overlap_EmpiricalCompounds:
                EIDs.append(EID_set.EID)
                cpds.append(EID_set.chosen_compounds)

            print(P.adjusted_p)

            names = [[self.pathway_analysis.mixedNetwork.model.dict_cpds_def.get(x, '') for x in y] for y in cpds]
            resultstr.append([str(x) for x in [P.name, P.overlap_size, P.EmpSize, P.adjusted_p]]
                             + [','.join(EIDs), ','.join(['/'.join(x) for x in cpds]),
                                '$'.join(['/'.join(x) for x in names])])

        self.pathway_enrich_df = self.pathway_enrich_df.from_records(resultstr)
        self.pathway_enrich_df.columns = self.pathway_enrich_df.iloc[0]
        self.pathway_enrich_df = self.pathway_enrich_df.iloc[1:]

        return self.pathway_enrich_df.copy()

        # try:
        # self.pathway_enrich_df.to_csv("pathway_enrichment_table" +
        # self.PathwayAnalysis.mixedNetwork.data.paradict['output'], sep = '\t')

        # except IOError:
        # print("specified path incorrect")

    def construct_mwas_plots(self):

        figsize = (15, 8)
        CutoffLine = -np.log10(self.PathwayAnalysis.mixedNetwork.data.paradict['cutoff'])

        sigList = [f for f in self.PathwayAnalysis.mixedNetwork.data.ListOfMassFeatures
                   if f.p_value < self.PathwayAnalysis.mixedNetwork.data.paradict['cutoff']]
        restList = [f for f in self.PathwayAnalysis.mixedNetwork.data.ListOfMassFeatures
                    if f.p_value >= self.PathwayAnalysis.mixedNetwork.data.paradict['cutoff']]

        Y_label = "-log10 p-value"
        Y_black = [-np.log10(f.p_value) for f in restList]
        Y_green = [-np.log10(f.p_value) for f in sigList]
        X_label = ["m/z", "Retention time"]
        X_black = [[f.mz for f in restList], [f.retention_time for f in restList]]
        X_green = [[f.mz for f in sigList], [f.retention_time for f in sigList]]
        X_max = [self.PathwayAnalysis.mixedNetwork.data.max_mz,
                 self.PathwayAnalysis.mixedNetwork.data.max_retention_time]

        fig, myaxes = plt.subplots(figsize=figsize, nrows=1, ncols=2)
        for ii in range(2):
            myaxes[ii].scatter(X_black[ii], Y_black, s=10, c='black', linewidths=6, alpha=0.8)
            myaxes[ii].scatter(X_green[ii], Y_green, s=10, c='green', linewidths=6, alpha=0.8)
            # lines
            myaxes[ii].plot([0, X_max[ii]], [CutoffLine, CutoffLine], 'g--')

            myaxes[ii].spines['right'].set_visible(True)
            myaxes[ii].spines['top'].set_visible(True)
            myaxes[ii].yaxis.set_ticks_position('both')
            myaxes[ii].xaxis.set_ticks_position('both')

            myaxes[ii].set_xlabel(X_label[ii])
            myaxes[ii].set_ylabel(Y_label)
            # rotate to avoid overlap xticklabels
            plt.setp(myaxes[ii].get_xticklabels(), rotation=30, horizontalalignment='right')

            myaxes[ii].xaxis.label.set_size(25)
            myaxes[ii].yaxis.label.set_size(25)

        #plt.tight_layout()
        plt.show()

        # return fig

    def pathway_significance_plot(self):

        self.PathwayAnalysis.permutation_record.sort()

        Y_data = [-np.log10(x) for x in self.PathwayAnalysis.permutation_record]
        fig = plt.figure(figsize=(5, 4))
        plt.plot(range(len(Y_data)), Y_data, 'b.')
        for P in self.PathwayAnalysis.resultListOfPathways[:10]:
            YY = -np.log10(P.p_EASE)
            plt.plot([0, 0.1 * len(Y_data)], [YY, YY], 'r--')

        plt.ylabel("-log10 (FET p-value)")
        plt.xlabel("Number of permutations")
        plt.title("Pathway significance")
        plt.tight_layout()

    def export_EmpiricalCompounds(self):

        s = "EID\tmassfeature_rows\tstr_row_ion\tcompounds\tcompound_names\n"
        for E in self.PathwayAnalysis.mixedNetwork.ListOfEmpiricalCompounds:
            names = [self.PathwayAnalysis.mixedNetwork.model.dict_cpds_def.get(x, '') for x in E.compounds]
            s += '\t'.join([E.EID, ';'.join(E.massfeature_rows), E.str_row_ion, ';'.join(E.compounds), '$'.join(names)]
                           ) + '\n'

        s = "input_row\tEID\tstr_row_ion\tcompounds\tcompound_names\tinput_row\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user\n"

        for row in self.PathwayAnalysis.mixedNetwork.mzrows:
            # not all input rows match to an empCpd
            try:
                for E in self.PathwayAnalysis.mixedNetwork.rowindex_to_EmpiricalCompounds[row]:
                    names = [self.PathwayAnalysis.mixedNetwork.model.dict_cpds_def.get(x, '') for x in E.compounds]
                    s += '\t'.join([row, E.EID, E.str_row_ion, ';'.join(E.compounds), '$'.join(names)]
                                   ) + '\t' + self.PathwayAnalysis.mixedNetwork.rowDict[row].make_str_output() + '\n'
            except KeyError:
                pass

        cols = s.split("\t")

        # print(cols)

        empc_df = pd.DataFrame(columns=['input_row', 'EID', 'str_row_ion',
                                        'compounds', 'compound_names', 'input_row', 'm/z', 'retention_time', 'p_value',
                                        'statistic', 'CompoundID_from_user'])

        print(empc_df.head())


class MummichogPathwayAnalysis(Method):

    def __init__(self, mummichog_data_source):
        # mummichog_data_source - Input User data object

        # logger.debug('Mummichog initialised')
        self.comparisons_list = []
        self.mummichog_datasource = mummichog_data_source

        for data_source in mummichog_data_source:
            self.theoretical_model = metabolicNetwork(metabolicModels['human_model_mfn'])
            self.mixed_network = DataMeetModel(self.theoretical_model, data_source)
            self.PA = PathwayAnalysis(self.mixed_network.model.metabolic_pathways, self.mixed_network)
            self.comparisons_list.append(self.PA)

        # print("result list output", self.PA.resultListOfPathways)

        # Modular Analysis and Activity network omitted change local run to a class method

    def get_results(self):
        #return a dataframe

        result_objects = []

        resultstr = [['pathway', 'comp1/comp2 pvalue',  'comp3/4 pvalue',
                      'overlap_EmpiricalCompounds (id)', 'overlap_features (id)', 'overlap_features (name)', ]]

        for comparison in self.comparisons_list:
            comparison.cpd_enrich_test()
            comparison.collect_hit_Trios()
            results = MummichogDataExporting(comparison)
            result_objects.append(results)



        return result_objects


    def get_mixed_network(self):
        return self.mixed_network

    def get_theoretical_model(self):
        return self.theoretical_model

        # perform the pathway analysis
