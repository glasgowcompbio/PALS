### 1. Introduction

**To access our interactive Web application *PALS Viewer*, please visit [https://pals.glasgowcompbio.org/app/](https://pals.glasgowcompbio.org/app/).**

Pathway analysis is an important task in understanding complex metabolomic data. Here we introduce **PALS (Pathway 
Activity Level Scoring)**, a complete tool that performs database queries of pathways, decomposes activity levels in pathways via [the PLAGE method](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-225), as well as presents the results in a user-friendly manner. The results are found to be more robust to noise and missing peaks compared to the alternatives (ORA, GSEA). This is particularly important for metabolomics peak data, where noise and missing peaks are prevalent.

![PALS](https://github.com/glasgowcompbio/PALS/raw/master/images/overall_schematic.png "PALS")

Additionally, the [decomposition approach](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-225) 
in PALS is amenable to the analysis of any group of metabolite sets, not just pathways. As demonstrated in PALS Viewer,
metabolite sets obtained from the grouping of metabolites according to their fragmentation spectra can also be analysed. 
This includes in particular *Molecular Families* from [GNPS](http://gnps.ucsd.edu/), as well as *Mass2Motifs* from [MS2LDA](http://ms2lda.org/). From [PALS Viewer](https://pals.glasgowcompbio.org/app/), you can also prioritise MF and Mass2Motifs from your GNPS analysis based on their activity levels. For more details including how to analyse user-defined metabolite sets from Jupyter notebooks, see Section 8.

### 2. Installation

For the latest development version, check out this repository using Git.
Otherwise PALS can also be installed via `pip install PALS-pathway`.

To use Reactome as pathway database, refer to the [setup guide](https://github.com/glasgowcompbio/PALS/blob/master/setup_guide.md).

### 3. Running PALS

To run PALS from the command-line, the script *pals/run_pals.py* is used. This script accepts a number of parameters, 
documented below (**bold** indicates required parameters). 

**Note: if you have installed PALS via pip, then *run_pals.py* (Unix-based systems) or
*run_pals.exe* (Windows) will also be added to your path during installation. 
It can be run directly by typing its name from the shell.**

#### Usage Parameters

```
PLAGE/ORA/GSEA
usage: run.py [-h] -x [--method] {PLAGE, ORA, GSEA} 
              
              
        --db, The pathway database to use. Valid choices are as follows. 
              - *PiMP_KEGG*: KEGG compound database exported from PiMP.
              - *COMPOUND*: Reactome compound database matching by KEGG ids.
              - *ChEBI*: Reactome compound database matching by ChEBI ids.
              - *UniProt*: Reactome protein database matching by UniProt ids.
              - *ENSEMBL*: Reactome gene database matching by ENSEMBL ids.        
              
              Note that *PiMP_KEGG*, *COMPOUND* and *ChEBI* are for metabolomics use, while *UniProt* and *ENSEMBL* are for 
              proteomics and transcriptomics use respectively. {PiMP_KEGG,COMPOUND,ChEBI,UniProt,ENSEMBL}
              
       --comparisons, Specifies the comparisons to make, e.g. `--comparisons beer1/beer2 beer3/beer4` to specify 
                      beer1 (case) vs beer2 (control), as well as beer3 (case) vs beer4 (control). COMPARISONS [COMPARISONS ...]
              
       --min_replace, The minimum intensity value for data imputation, e.g. `--min_replace 5000`. Defaults to 5000.      
              
       --species, Species name for Reactome pathway query, e.g. `--species "Homo sapiens"`. Defaults to Homo Sapiens. 
                  {Arabidopsis thaliana,Bos taurus,Caenorhabditis elegans,
                   Canis lupus familiaris,Danio rerio,Dictyostelium discoideum,
                   Drosophila melanogaster,Gallus gallus,Homo sapiens,Mus musculus,
                   Oryza sativa,Rattus norvegicus,Saccharomyces cerevisiae,Sus scrofa}
            
              [--use_all_reactome_pathways] [--connect_to_reactome_server]
              intensity_csv annotation_csv output_file

MUMMICHOG
usage: run_pals.py [-h] -x [--method] {MUMMICHOG}
       
       Basic usage example:
       python run_pals.py -f mydata.txt -o myoutput
       
       -f, --infile: single file as input, 
              containing all features with tab-delimited columns
              m/z, retention time, p-value, statistic score

        -n, --network: network model to use (default human_mfn; models being ported to version 2), 
              [human_mfn, worm]

        -o, --output: output file identification string (default 'mcgresult')
        -k, --workdir: directory for all data files.
              Default is current directory.

        -m, --mode: analytical mode of mass spec, [positive, negative, pos_defult].
              Default is pos_defult, a short version of positive.
        -u, --instrument: Any integer, treated as ppm of instrument accuracy. Default is 10. 

        -p, --permutation: number of permutation to estimate null distributions.
              Default is 100.
        -z,   --force_primary_ion: one of primary ions, 
              ['M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]'],  
              must be present for a predicted metabolite, [True, False].
              Default is True.

        -c, --cutoff: optional cutoff p-value in user supplied statistics,
              used to select significant list of features. 
        -d, --modeling: modeling permutation data, [no, gamma].
              Default is no.
              
```
- **method**: Pathway ranking method to use, e.g. PLAGE, ORA, GSEA, MUMMICHOG.
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

#### Command-line Examples

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

Users provide two input files to PALS. The first is a matrix is of individual peak intensities (rows are peak features with 
column one containing the peak id, further columns representing individual samples). If using the CSV file, the second 
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

Data imputation is performed to the intensity matrix when it is loaded: if all of the samples in a single experimental factor 
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

ORA and GSEA are used as comparisons and are both included in PALS, and can be used for benchmarking and analysis. For more details, please refer to our paper.

We have made available a topological analysis method; mummichog in PALS. This requires a single tab-separated table containing mass-to-charge (m/Z) and retention time data in its usage.

### 6. Web Interface

PALS Viewer is a Web interface on top of the [Streamlit](https://www.streamlit.io/) framework. It can be used to run PALS, analyse 
pathway ranking results as well as inspect significantly changing pathways. To run it locally use following command:
```
$ streamlit run pals/run_gui.py
```
Alternatively an instance of PALS Viewer can also be found on our server.

### 7. Importing PALS for use in other applications

PALS can be imported as a Python library and incorporated into your own Python application. This is illustrated in the 
following code snippet:
```
from pals.ORA import ORA
from pals.PLAGE import PLAGE
from pals.GSEA import GSEA
from pals.common import *
from pals.feature_extraction import DataSource

# TODO: correctly initialise the following data structures for your data.
# Refer to the paragraphs below.
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

method = PLAGE(ds)
# method = ORA(ds) 
# method = GSEA(ds)

df = method.get_pathway_df()    
```

When PALS is used programatically, pandas dataframes storing relevant data can be passed directly. Experimental design 
data can be passed directly as a dictionary structure in the programmatic use. In the example above, `int_df` is the 
intensity data frame containing peak intensity information described in the File Format section above (with the second 
line of grouping information omitted). Similarly `annot_df` is the annotation data frame containing peak annotations 
also described above. 

Example `int_df` dataframe:
![int_df](https://github.com/glasgowcompbio/PALS/raw/master/images/int_df.png "int_df")

Example `annot_df` dataframe:
<img src="https://github.com/glasgowcompbio/PALS/raw/master/images/annot_df.png" width="200">

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
Many example notebooks are provided in the [notebooks](https://github.com/glasgowcompbio/PALS/tree/master/notebooks) folders.
If you have any question in using PALS in your application, please raise an issue or simply email us.'

### 8. Analysis of User-defined Metabolite Sets

PALS can be used to analyse any user-defined group of metabolites. In the paper, we demonstrated this through an example
by analysing potentially unknown metabolites that have been grouped into Molecular Families (by GNPS) and into
Mass2Motifs (by MS2LDA). [The following notebook](https://github.com/glasgowcompbio/PALS/blob/master/notebooks/GNPS_analysis.ipynb) demonstrates how the analysis was done.
Additionally Molecular Families and Mass2Motifs analysis have also been included as an option in PALS Viewer.

### 9. Publication

If you are using PALS, please cite the following publication.

Mcluskey, K., Wandy, J., Vincent, I., Hooft, J. J. J. Van Der, Rogers, S., Burgess, K. & Daly, R. (2021). *Ranking Metabolite Sets by Their Activity Levels.* Metabolites. 11(2):103 https://doi.org/10.3390/metabo11020103
