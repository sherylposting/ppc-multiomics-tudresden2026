# generates tiny_cpg, first 20k bp of each sample to do quick testing with pipeline

for i in /projects/p_dna15/data/EM_seq_files/cytosine_reports/*.gz; do 
    filename=$(basename "$i")
    zcat $i | head -n 20000 | gzip > "/projects/p_dna15/data/EM_seq_files/cytosine_reports/tiny_cpg/tiny_$filename"; 
done