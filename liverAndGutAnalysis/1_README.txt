####IMPORTANT NOTE!!!!!! Some paths in these scripts will point to /home/jovyan/Dediff (local jupyterhub environment). before re-running, check that these are changed to  /lustre/scratch117/casm/team274/gk14/Dediff. The two folders have been R-synced, so as long as the path before Dediff is changed, it should work fine.####

1. 1_stemness_ee_cor.R - script to calculate/plot correlation between early embryo and stemness signal, and also the change in early embryo signal when gaatrulation reference is added (Figs 4.1 and 4.2 in my thesis). Matt calculated the values used in the script. 
2. Scripts 1_gut_cellsig_file_setup.R and 4_liver_cellsig_file_setup.R are scripts to set up all the things I initially needed for tissue-specific cell signal analysis - processing reference data, reannotating or merging annotations where needed, etc. Also saving bulk RNA-seq files in the format needed for cell signal analysis. reformat_blastoid_samples.R is where I save blastoid samples (i.e. bulk samples that should match early embryo) into the desired format. Single-cell outputs of these files are used in the make_ee_gast_cellsig_refs.R, where I make the single-cell sumamry files for early embryo, early embryo + gastrulation, early embryo + gastrulation + tissue-specific references. Unless something needs reprocessing/reannotating, you don't really need to run the file_setup scripts anymore as all relevant outputs are in /lustre/scratch117/casm/team274/gk14/Dediff/. 
3. make_ee_gast_cellsig_refs.R - as previously mentioned, this script makes "final" reference sumamry files used in cell signal analysis, including processing gastrulation data which is smartseq2.  Takes some of the outputs of cellsig_file_setup ( see above ) scripts as inputs. 
4. run_cellsig_all.sh contains commands to run all cell signal analysis needed for this project. e.g. each bulk type (liver, gut) against each reference (ee, ee + gast, etc etc). Don't forget to chage all instances of /home/jovyan/Dediff...
5. Scripts make_gut_bulk_metadata.R and make_liver_bulk_metadata_file.R are pretty self-explanatory. They make a "metadata" file which contains metadata file for all paths in the -b cell signal analysis input file. These metadata files are saved, so no need to re-run these scripts. If you end up updating your bulk files, you'll want to re-make the metada yourself anyway as this is a fairly manual process...
6. Scripts 3_gut_single_cell.R and 6_liver_single_cell.R contain the setup and LR analysis of tumor single-cell against complete refs (i.e. ee+gast+tissue-specific). 






