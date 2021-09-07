from pals import feature_extraction
from pals.common import load_data, SMALL, DATABASE_PIMP_KEGG
import pandas as pd
import numpy as np

from PALS.pals.RunMummichog import MummichogPathwayAnalysis
from PALS.pals.feature_extraction import DataSource

intdf, annodf, groups = load_data("int_df.csv", "annotation_df.csv") #path to beer data

mz_df = pd.DataFrame(columns = ['m/z', 'retention_time'])

mz_df['m/z'] = [np.random.randint(1, 150) for i in range(7375)]
mz_df['retention_time'] = [np.random.randint(1, 300) for i in range(7375)]


comparisons = [
    ('beer1', 'beer2'),
    ('beer3', 'beer4')
]

experimental_design = {
    'groups': groups,
    'comparisons': []
}
for case, control in comparisons:
    experimental_design['comparisons'].append({
        'case': case,
        'control': control,
        'name': '%s/%s' % (case, control)
    })


ds = DataSource(intdf, annodf, experimental_design, DATABASE_PIMP_KEGG, min_replace=SMALL, mz_rt_df= mz_df, comparisons= comparisons)

mummichog = ds.create_mummichog_ds()

pathway_analysis = MummichogPathwayAnalysis(mummichog)
data_objects = pathway_analysis.get_results()

pathway_df = data_objects[1].construct_pathway_enrich_df()

pathway_df.to_csv("new2_df.csv")