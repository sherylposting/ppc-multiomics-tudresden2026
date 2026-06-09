- use `ls -a` to display hidden files (such as .git, .gitignore)
- git is ignoring raw data files (data/)
- code assumes your working directory is set to project root

### Directory structure
/projects/p_dna15 \\
|  data : contains all raw data \\
|  |  EM_seq_files
|  |  |  corr.txt
|  |  |  bismark .cov.gz files...
|  |  |  BFX2892-EM-Seq_multiqc_report
|  |  |  L188015_WT1_R1_bismark_bt2_pe.deduplicated_splitting_report
|  |  |  L188015_WT1_R1_bismark_bt2_PE_report
|  code
|  |  preliminary_v1.0.R
|  methylDB : files needed for tabix DB
|  .git
|  .gitignore
|  README.md
