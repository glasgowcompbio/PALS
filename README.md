# PALS
Pathway-level Analysis of Metabolites Expression Data

`pip install PALS-pathway`

Command-line Usage
------------------

Runs PALS locally (not connecting to Reactome)

$ `python run_pals.py ../notebooks/test_data/beer/int_df.csv ../notebooks/test_data/beer/annotation_df.csv ../notebooks/test_data/beer/test_output.csv --db PiMP_KEGG --comparisons beer1/beer2 beer3/beer4 --min_replace 5000 --plage_weight 5 --hg_weight 1`

Runs PALS using KEGG Reactome on metabolic pathways only from Homo sapiens

$ `python run_pals.py ../notebooks/test_data/beer/int_df.csv ../notebooks/test_data/beer/annotation_df.csv ../notebooks/test_data/beer/test_output.csv --db COMPOUND --comparisons beer1/beer2 beer3/beer4 --min_replace 5000 --plage_weight 5 --hg_weight 1 --species "Homo sapiens"`

Runs PALS using KEGG Reactome on all pathways from Homo sapiens

$ `python run_pals.py ../notebooks/test_data/beer/int_df.csv ../notebooks/test_data/beer/annotation_df.csv ../notebooks/test_data/beer/test_output.csv --db COMPOUND --comparisons beer1/beer2 beer3/beer4 --min_replace 5000 --plage_weight 1 --hg_weight 1 --species "Homo sapiens" --use_all_reactome_pathways --connect_to_reactome_server`
