# PALS

### Motivation

Understanding changing pathways can be incredibly useful in the interpretation and understanding of complex datasets from metabolomics experiments. While there is an abundance of tools for enrichment analysis of genomics data, those suitable for metabolomic pathway analysis are still relatively scarce.  Many of the available metabolomic tools perform simple over-representation analysis and do not fully utilise valuable mass-spectrometry peak intensity data. In addition, it is often difficult to integrate these tools into a Python-based workflow and incorporate complex experimental designs.

### Results

Here we present **PALS (Pathway Activity Level Scoring)**, a Python package to perform the ranking of significantly-changing metabolite pathways in different experimental conditions. PALS achieves this through the decomposition of pathway activity levels calculated from peak intensities. PALS can operate either in an offline mode by using previously downloaded KEGG and Reactome databases included in the package, or in an online mode by utilising a local Reactome service to retrieve the most up-to-date pathways. Results from running a PALS analysis can be exported as a comma-separated file or viewed interactively in an online Web platform. A comparison of PALS with two other commonly used methods for metabolomic pathway analysis (ORA and GSEA) is also given, and reveals that PALS is more robust to missing peaks and noisy data than the alternatives. Additionally, PALS is used to analyse pathways from a study of Human African Trypanosomiasis and the results reported.

### Installation

For the latest bleeding-edge version, check out this repository using Git.
Otherwise PALS can also be installed via `pip install PALS-pathway`.

To use Reactome as pathway database, refer to the [setup guide](setup_guide.md).

### Command-line Usage

The most basic usage of PALS is to run it in offline-mode, which uses the downloaded KEGG database for pathways. Here we run PALS on a test beer data.

***TODO: replace this with the HAT data in the manuscript***

$ `python pals/run.py PALS notebooks/test_data/beer/int_df.csv notebooks/test_data/beer/annotation_df.csv test_output.csv --db PiMP_KEGG --comparisons beer1/beer2 beer3/beer4`

Downloaded Reactome pathways is also provided in PALS for the [most common species](https://github.com/glasgowcompbio/PALS/tree/master/pals/data/reactome/metabolic_pathways/COMPOUND). Note that only metabolic pathways are available in this mode. Below shows an example run using Human pathways.
 
$ `python pals/run.py notebooks/test_data/beer/int_df.csv notebooks/test_data/beer/annotation_df.csv test_output.csv --db COMPOUND --comparisons beer1/beer2 beer3/beer4 --min_replace 5000 --species "Homo sapiens"`

Finally in online mode, PALS can connect to a Reactome database instance to retrieve the most updated pathways. The flag `--use_all_reactome_pathways` specifies that all pathways should be used (not just metabolic pathways), while the flag `--connect_to_reactome_server` defines that online mode should be used.

$ `python pals/run.py notebooks/test_data/beer/int_df.csv notebooks/test_data/beer/annotation_df.csv test_output.csv --db COMPOUND --comparisons beer1/beer2 beer3/beer4 --min_replace 5000 --species "Homo sapiens"  --use_all_reactome_pathways --connect_to_reactome_server`

### Example Results

### Other Pathway Analysis Methods

ORA, GSEA ..