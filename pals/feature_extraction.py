import os
import xmltodict

from loguru import logger


class DataSource(object):

    def __init__(self, int_df):
        self.analysis = analysis
        self.project = Analysis.get_project(analysis)
        self.dataset = analysis.dataset_set.all()[0]
        self.comparisons = list(Comparison.objects.filter(experiment=analysis.experiment))
        self.path = os.getcwd()

        # a dataframe of peak intensities, where rows = ms1_peak_id and columns = sample_name
        self.int_df = self.construct_peak_int_df()

        # Dictionaries required throughout the class
        self.pathway_cmpd_dict = self.get_pw_cmpd_dict()
        self.cmpd_formula_dict = self.get_cmpd_name_formula_dict()
        self.cmpd_id_formula_dict = self.get_cmpd_id_formula_dict()
        self.cmpd_id_pw_dict = self.get_cmpid_pathway_dict()

        # For use in the hypergeometric test - the number of unique formulas in kegg and in pathways
        # and the number of unique formulas in the ds and in pathways
        self.kegg_pw_formulas, self.unique_ds_pw_fs = self.get_unique_pw_f()
        self.groups = self.get_groups()

    """A method to construct a peak_intestity dataframe with rows = peak_ids and columns = sample names
    Adapted from experiments.views.get_peak_table
    :returns: DF with index of peak Ids, columns of sample names and values of peak intensites.
    """

    def construct_peak_int_df(self):
        logger.info("Constructing the peak intensity DF")
        recs = PeakDTSample.objects.filter(peak__dataset=self.dataset).values_list('peak', 'sample', 'intensity')
        peaks, samples, intensities = zip(*recs)

        dpeaks = sorted(list(set(peaks)))
        dsamples = sorted(list(set(samples)))

        peakdict = dict(zip(dpeaks, range(len(dpeaks))))
        sampledict = dict(zip(dsamples, range(len(dsamples))))

        rows = [peakdict[_] for _ in peaks]
        cols = [sampledict[_] for _ in samples]

        data = coo_matrix((intensities, (rows, cols)))
        sample_names = Sample.objects.filter(id__in=dsamples).order_by('id').values_list('name', flat=True)
        int_df = pd.DataFrame(data.todense(), index=dpeaks, columns=sample_names)
        int_df.index.name = "ms1_peak_id"
        int_df.columns.name = "sample_name"

        return int_df

    """Returns: dictionary pathway (mapid) : [compounds in that pathway (by name)]
    """

    def get_pw_cmpd_dict(self):

        robjects.r['load'](self.path + '/PiMP/data/pathways2Compounds.RData')
        a = robjects.r['pathways2Compounds']
        pw_cmpd_dict = OrderedDict(zip(a.names, map(list, list(a))))

        pathway_cmpd_dict = OrderedDict()
        for pathway, cmpd_list in pw_cmpd_dict.iteritems():
            pw = pathway.replace("path:", "")
            pathway_cmpd_dict[pw] = cmpd_list

        return pathway_cmpd_dict

    """ Returns: a dictionary with the compound_name:formula
    """

    def get_cmpd_name_formula_dict(self):

        return self._produce_kegg_dict('name')

    """ Returns: a dictionary with the compound_id:formula
    """

    def get_cmpd_id_formula_dict(self):

        return self._produce_kegg_dict('id')


    """ Returns: a dictionary with the compound_id:[pathways]
    """

    def get_cmpid_pathway_dict(self):

        robjects.r['load'](self.path + '/PiMP/data/compounds2Pathways.RData')
        a = robjects.r['compounds2Pathways']
        cmpd_pw_dict = OrderedDict(zip(a.names, map(list, list(a))))

        return cmpd_pw_dict

    """
    Method to return the unique pathway formulas associated with an analysis
    """

    def get_unique_pw_f(self):
        # set of compounds by id in Kegg
        kegg_cmpd_ids = set(self.cmpd_id_formula_dict.keys())
        # set of compounds by id found in the pathways
        pathway_cmpd_ids = set(self.cmpd_id_pw_dict.keys())
        # These are the pathway compound (by id) that are also present in Kegg
        pw_cmpd_ids = kegg_cmpd_ids.intersection(pathway_cmpd_ids)

        kegg_pathway_formulas = {}
        for c in pw_cmpd_ids:
            kegg_pathway_formulas[c] = self.cmpd_id_formula_dict[c]
        kegg_pw_formulas = set(kegg_pathway_formulas.values())

        # ds_formulas that are M+H or M-H - probably should use
        ds_formulas = Compound.objects.filter(Q(identified=True) | (Q(identified=False) & Q(adduct__in=["M+H", "M-H"])),
                                              repositorycompound__db_name='kegg',
                                              peak__dataset=self.dataset).values_list(
            "formula", flat=True).distinct()

        unique_ds_pw_fs = kegg_pw_formulas.intersection(set(ds_formulas))
        return len(kegg_pw_formulas), len(unique_ds_pw_fs)

    """ Method to return a dictionary of the groups (factors) and sample names
    (adapted from Joe's code)
        Returns:groups_dict: A dictionary of the group (key) and samples in that group (values)
    """

    def get_groups(self, int_df):

        m = Rpy2PipelineMetadata(self.analysis, self.project)
        groups = m.get_groups(strip_extension=False)
        groups_df, _ = convert_to_dataframe(groups)
        # Filter out those samples that are not included in the analysis
        groups_df = groups_df.loc[groups_df['sample'].isin(int_df.columns)]

        groups_dict = defaultdict(list)
        for idx, row in groups_df.iterrows():
            sample_name = row.values[0]
            factors = tuple(row.values[1:-1])
            groups_dict[factors].append(sample_name)

        return groups_dict, groups_df

    """ Param: a comparison object
        Returns: a list of lists containing the file names of the samples
        for each attribute in the comparison
    """

    def get_comparison_samples(self, comparison):
        attribute = Attribute.objects.filter(comparison=comparison)
        samples = []
        for a in attribute:
            samples.append(Sample.objects.filter(attribute=a).values_list('name', flat=True))
        return samples

    """ Method to get the comparison attributes as a string <Condition>/<Control>
        Param: a comparison object
        Returns: a string containing the name of the attributes compared in the comparison
    """

    def get_comparison_names(self, comparison):
        att_group = []
        attributeC = AttributeComparison.objects.filter(comparison=comparison).distinct().order_by('id')
        for a in attributeC:
            att_group.append((a.attribute.name, a.group,))

        sort_groups = sorted(att_group, key=lambda x: x[1])
        comp_names = sort_groups[1][0] + "/" + sort_groups[0][0]
        return comp_names

    """ Param: the mapid of a pathway (e.g. map00010)
        Returns: the number of unique formula identifiable pathway
    """

    def get_pw_unique_F(self, mapid):
        self.pathway_cmpd_dict
        self.cmpd_formula_dict
        compounds = self.pathway_cmpd_dict[mapid]
        formulae = set()
        for c in compounds:
            try:
                formula = self.cmpd_formula_dict[c]
                formulae.add(formula)
            except KeyError:
                continue
        return len(formulae)

    """ Param: list of compound ids
        Returns: the number of unique formula in the list
    """

    # New method to tie in with get_pathway_compounds from experiment.views which returns compound ids
    def get_unique_F_by_id(self, compound_ids):
        c_ids = compound_ids.split(',')
        numberFormulae = Compound.objects.filter(id__in=c_ids).values_list('formula', flat=True).distinct().count()
        return numberFormulae

    """
    Takes in a dataset and p-value adjustment and returns a list of lists containing all the pathways and associated compounds.
    [[pathway primary key, comma separated string of compound ids in that pathway, pathway name], ...]
    Moved from views so not to have circular calls but also think it makes sense in here (KMcL)
    """

    def get_ds_pw_compounds(self, p_val_adjustment='none'):
        if p_val_adjustment == 'corrected':
            peaks = Peak.objects.filter(dataset=self.dataset, peakdtcomparison__adjPvalue__lt=0.05).distinct()
        elif p_val_adjustment == 'uncorrected':
            peaks = Peak.objects.filter(dataset=self.dataset, peakdtcomparison__pValue__lt=0.05).distinct()
        elif p_val_adjustment == 'none':
            peaks = Peak.objects.filter(dataset=self.dataset).distinct()

        # A compound can be in more than one pathway
        compounds = Compound.objects.filter(Q(identified=True) | (Q(identified=False) & Q(adduct__in=["M+H", "M-H"])),
                                            repositorycompound__db_name='kegg', peak__in=peaks)

        pathways = Pathway.objects.filter(datasourcesuperpathway__compoundpathway__compound__in=compounds).values_list(
            'id', 'datasourcesuperpathway__compoundpathway__compound__id', 'name')

        p_dict = OrderedDict()
        for p in pathways:
            if p[2] in p_dict:
                p_dict[p[2]]['compound_ids'].append(p[1])
            else:
                p_dict[p[2]] = {'compound_ids': [p[1]],
                                'id': p[0]}

        for k in p_dict.keys():
            p_dict[k]['compound_ids'] = ','.join([str(c_id) for c_id in p_dict[k]['compound_ids']])

        pathway_info = [[p_info['id'], p_info['compound_ids'], p_name] for p_name, p_info in p_dict.items()]
        return pathway_info

    ####################################################################################################################
    # private methods
    ####################################################################################################################

    """Param:
        Parameter to define what keys we want from the kegg_cmpd_dict
    Returns:
        a dictionary from kegg.xml depending on the parameter that is passed
    """

    def _produce_kegg_dict(self, param):
        kegg_location = self.path + '/PiMP/inst/dbs/kegg.xml'
        with open(kegg_location) as kegg_cmpd_file:
            cmpd_dict = xmltodict.parse(kegg_cmpd_file.read())

        kegg_dict = {}
        for compound in cmpd_dict['compounds']['compound']:
            kegg_dict[compound[param]] = compound['formula']
        return kegg_dict

