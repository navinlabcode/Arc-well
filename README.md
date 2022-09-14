# Arc-well
This repository contains the scripts, meta data metrics files used in the manuscript: Archival single cell sequencing reveals persistent subclones over years to decades of DCIS progression.


##### scripts 
includes all of the codes that could reproduce figures of Arc-well paper.
- _snakemake_: includes the workflow of downsampling coverage calculations
- _CNA_pipeline_: includes the workflow of initial QC and variable binning, prepares the inputs for downstream analysis

##### pre_load_data 
files that required during excuting the R scripts under scripts folder
##### metrics  
includes filtering files, QC metrics, clinical, computational meta data from cell lines, FFPE samples.

#### _Docker image_
------------
In case you have problem in installing certain version of dependencies, you can run scripts without installing dependencies by pull the docker image here.
#### _Data source_
------------
The original sequencing data from this study has been deposited to the Sequence Read Archive (SRA): PRJNA799605.

#### _Contact_
------------
For more detailed information, please email nnavin@mdanderson.org 

