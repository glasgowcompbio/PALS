# PALS

<iframe width="560" height="315" src="https://www.youtube.com/embed/XP4wsoInh4E" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

### 1. Abstract

Understanding changing pathways can be incredibly useful in the interpretation and understanding of complex datasets from 
metabolomics experiments. While there is an abundance of tools for enrichment analysis of genomics data, those suitable 
for metabolomic pathway analysis are still relatively scarce.  Many of the available metabolomic tools perform simple 
over-representation analysis and do not fully utilise valuable mass-spectrometry peak intensity data. In addition, it 
is often difficult to integrate these tools into a Python-based workflow and incorporate complex experimental designs.

Here we present **PALS (Pathway Activity Level Scoring)**, a Python library to perform the ranking of significantly-changing 
metabolite pathways in different experimental conditions. PALS achieves this through the decomposition of pathway activity 
levels calculated from peak intensities. PALS can operate either in an offline mode by using previously downloaded KEGG 
and Reactome databases, or in an online mode by utilising a local Reactome service to retrieve the most up-to-date pathways. 

![PALS](images/overall_schematic.png?raw=true "PALS")

Results from running a PALS analysis can be exported as a comma-separated file or viewed interactively in an online Web platform. 
A comparison of PALS with two other commonly used methods for metabolomic pathway analysis (ORA and GSEA) is also given, 
and reveals that PALS is more robust to missing peaks and noisy data than the alternatives. Additionally, PALS is used 
to analyse pathways from a study of Human African Trypanosomiasis and the results reported.

### 2. Installation

For the latest bleeding-edge version, check out this repository using Git.
Otherwise PALS can also be installed via `pip install PALS-pathway`.

To use Reactome as pathway database, refer to the [setup guide](setup_guide.md).

### 3. Command-line Usage

To run PALS from the command-line, the script *pals/run.py* is used. This script accepts a number of parameters, 
documented here (**bold** indicates required parameters):
```
usage: run.py [-h] --db {PiMP_KEGG,COMPOUND,ChEBI,UniProt,ENSEMBL}
              --comparisons COMPARISONS [COMPARISONS ...]
              [--min_replace MIN_REPLACE]
              [--species {Arabidopsis thaliana,Bos taurus,Caenorhabditis elegans,
              Canis lupus familiaris,Danio rerio,Dictyostelium discoideum,
              Drosophila melanogaster,Gallus gallus,Homo sapiens,Mus musculus,
              Oryza sativa,Rattus norvegicus,Saccharomyces cerevisiae,Sus scrofa}]
              [--use_all_reactome_pathways] [--connect_to_reactome_server]
              {PLAGE,ORA} intensity_csv annotation_csv output_file
```
- **method**: Pathway ranking method to use, e.g. PLAGE or ORA.
- **intensity_csv**: Input intensity CSV file
- **annotation\_csv**: Input annotation CSV file
- **output_file**: Output pathway ranking file.
- **--db**: The pathway database to use. Valid choices are as follows. 
    - *PiMP_KEGG*: KEGG compound database exported from PiMP.
    - *COMPOUND*: Reactome compound database matching by KEGG ids.
    - *ChEBI*: Reactome compound database matching by ChEBI ids.
    - *UniProt*: Reactome protein database matching by UniProt ids.
    - *ENSEMBL*: Reactome gene database matching by ENSEMBL ids.        
    Note that *PiMP_KEGG*, *COMPOUND* and *ChEBI* are for metabolomics use, while *UniProt* and *ENSEMBL* are for 
    proteomics and transcriptomics use respectively.
- **--comparisons**: Specifies the comparisons to make, e.g. `--comparisons beer1/beer2 beer3/beer4` to specify 
beer1 (case) vs beer2 (control), as well as beer3 (case) vs beer4 (control).
- --min_replace: The minimum intensity value for data imputation, e.g. `--min_replace 5000`. Defaults to 5000.
- --species: Species name for Reactome pathway query, e.g. `--species "Homo sapiens"`. Defaults to Homo Sapiens.
- --use_all_reactome_pathways: Whether to use all pathways for Reactome pathway query. If this option is not used, 
only metabolic pathways will be queried.
- --connect_to_reactome_server: Whether to connect to an instance of Neo4j server hosting Reactome database 
(online mode). If not specified, then offline mode (using a downloaded copy of selected Reactome pathways) will be used. 

The most basic usage of PALS is to run it in offline-mode using PLAGE as the decomposition method. This uses the 
downloaded KEGG database for pathways. Here we run PALS on the example HAT data used in the manuscript

$ `python pals/run.py PLAGE notebooks/test_data/HAT/int_df.csv notebooks/test_data/HAT/annotation_df.csv test_output.csv --db PiMP_KEGG --comparisons Stage_1/Control Stage_2/Control`

Downloaded Reactome pathways is also provided in PALS for the 
[most common species](https://github.com/glasgowcompbio/PALS/tree/master/pals/data/reactome/metabolic_pathways/COMPOUND). 
Note that only metabolic pathways are available in this mode. Below shows an example run using Human pathways.
 
$ `python pals/run.py PLAGE notebooks/test_data/HAT/int_df.csv notebooks/test_data/HAT/annotation_df.csv test_output.csv --db COMPOUND --comparisons Stage_1/Control Stage_2/Control --min_replace 5000 --species "Homo sapiens"`

Finally in online mode, PALS can connect to a Reactome database instance to retrieve the most updated pathways. 
The flag `--use_all_reactome_pathways` specifies that all pathways should be used (not just metabolic pathways), while the flag `--connect_to_reactome_server` defines that online mode should be used.

$ `python pals/run.py PLAGE notebooks/test_data/HAT/int_df.csv notebooks/test_data/HAT/annotation_df.csv test_output.csv --db COMPOUND --comparisons Stage_1/Control Stage_2/Control --min_replace 5000 --species "Homo sapiens"   --use_all_reactome_pathways --connect_to_reactome_server`

The following example output is produced:

![output](images/output.png?raw=true "Output")

Pathways are identified by their id and can be sorted by the `p-value` columns. The column `unq_pw_F` lists the unique 
formulae found in that pathway, `tot_ds_F` lists the formula hits found in the dataset, and `F_coverage` is the proportion
of `tot_ds_F` to `unq_pw_F`.

### 4. File format and data imputation

The first input users have to provide is a matrix is of individual peak intensities where rows are peak features with 
column one containing the peak id and further columns representing individual samples. If using the CSV file, the second 
line can be used to indicate which groups this sample belongs to. For example, the intensity matrix takes the form of:
```
row_id,Beer_1_full1.mzXML,Beer_1_full2.mzXML,Beer_1_full3.mzXML,
    Beer_2_full1.mzXML,Beer_2_full2.mzXML,Beer_2_full3.mzXML,
    Beer_3_full1.mzXML,Beer_3_full2.mzXML,Beer_3_full3.mzXML,
    Beer_4_full1.mzXML,Beer_4_full2.mzXML,Beer_4_full3.mzXML
group,beer1,beer1,beer1,beer2,beer2,beer2,beer3,beer3,beer3,beer4,beer4,beer4
3033929,2235291136,2000478208,2170697216,2242759936,2279881984,1959479680,
    2079356160,2110473216,2243652608,1817064704,1746442752,1779827200
3033930,44334908,42873872,48948532,47604480,42172796,39084524,
    38257776,37701920,40871888,33304766,31536296,31024098
```

Data imputation is performed as such to the intensity matrix: if all of the samples in a single experimental factor 
have intensities of zero these are replaced by the minimum intensity value (which can be set by the user); and if 
only some of the sample values in a factor are zero then these are replaced by the mean value of the non-zero samples 
in that factor. The data is subsequently transformed to log-2 base and standardised using the preprocessing module in 
Scipy such that the intensity matrix has a zero mean and unit variance across the samples.

In addition, users also provide a list of compound annotations assigned to peak features (peaks that do not have 
annotations will not be used for pathway analysis). As a result of the uncertainty in peak identification, multiple 
peak IDs may be mapped to multiple compound IDs and \textit{vice versa}. As such, annotations are provided as 
another matrix having two columns. The first column (or DataFrame index) is the peak ID while the second column is 
the assigned metabolite annotation as either KEGG or ChEBI database IDs.
```
row_id,entity_id
3033929,C00148
3033930,C06326
3033931,C00183
3033931,C00719
3033931,C00583
```

### 5. Other Pathway Analysis Methods

ORA and GSEA are both included in PALS, and can be used for comparison. For more details, please refer to our paper.

### 6. Graphical User Interface

An experimental user interface on top of the [Streamlit](https://www.streamlit.io/) framework is also available and can 
be used as an alternative to the command-line interface.
![GUI](images/pals_viewer.png?raw=true "PALS Viewer")

PALS Viewer is a user-friendly graphical user interface to run PALS and analyse pathway ranking results as well as 
inspect significantly changing pathways. It can be run using the following command:

To run it, type:
```
$ streamlit run pals/run_gui.py
```

An online instance of PALS Viewer can also be accessed from ... **(TODO)**

### 7. Usage in other applications

PALS can be imported as a Python library and incorporated into your own Python application. This is illustrated in the 
following code snippet:
```
from pals.ORA import ORA
from pals.PALS import PALS
from pals.common import *
from pals.feature_extraction import DataSource

# TODO: correctly initialise the following data structures for your data.
int_df = pd.DataFrame()
annotation_df = pd.DataFrame()
experimental_design = {}

# Using Reactome pathways matching by KEGG ID
database_name = 'COMPOUND'

# If true, we limit to metabolic pathways only. Otherwise all pathways will be queried.
reactome_metabolic_pathway_only = True

# If true, we use online mode that queries Reactome on a local Neo4j server. 
# Otherwise offline mode will be used (using downloaded database files).
reactome_query = True

# Minimum intensity value for data imputation
min_replace = 5000

ds = DataSource(int_df, annotation_df, experimental_design, database_name,
                reactome_species=reactome_species,
                reactome_metabolic_pathway_only=reactome_metabolic_pathway_only,
                reactome_query=reactome_query, min_replace=min_replace)

method = PALS(ds)
# method = ORA(ds) 

df = method.get_pathway_df()    
```

When PALS is used programatically, pandas dataframes storing relevant data can be passed directly. Experimental design 
data can be passed directly as a dictionary structure in the programmatic use. In the example above, `int_df` is the 
intensity data frame containing peak intensity information described in the File Format section above (with the second 
line of grouping information omitted). Similarly `annot_df` is the annotation data frame containing peak annotations 
also described above. 

The experimental design data in `experimental_design` contains information on **groups**, which relates all samples in 
a particular experimental factor together as well as **comparisons**, which describes the desired comparisons for the 
PALS analysis in terms of a case and a control. An example of this can be found below:
```
experimental_design = {
    'comparisons': [
        {'case': 'beer1', 'control': 'beer2', 'name':  'beer1/beer2'},
        {'case': 'beer3', 'control': 'beer4', 'name': 'beer3/beer4'}
    ],
    'groups': {
        'beer1': [
            'Beer_1_full2.mzXML', 
            'Beer_1_full1.mzXML',
            'Beer_1_full3.mzXML'
        ],
        'beer2': [
            'Beer_2_full3.mzXML', 
            'Beer_2_full1.mzXML', 
            'Beer_2_full2.mzXML'
        ],
        'beer3': [
            'Beer_3_full3.mzXML', 
            'Beer_3_full2.mzXML', 
            'Beer_3_full1.mzXML'
        ],
        'beer4': [
            'Beer_4_full3.mzXML',
            'Beer_4_full2.mzXML',
            'Beer_4_full1.mzXML'
        ],
    }
}    
```