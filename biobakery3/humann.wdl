workflow humann {

	call run_humann

}

task run_humann {
	File fastq_R1
	File fastq_R2
	File fastq_S1
	File fastq_S2
	File metaphlan_profile_tsv
	String sample_id
	Float input_file_size = size(fastq_R1, "G") + size(fastq_R2, "G") + size(fastq_S1, "G") + size(fastq_S2, "G")
	
	command {    	

    	cat ${fastq_R1} ${fastq_R2} ${fastq_S1} ${fastq_S2} > input.fastq.gz
    	humann -i input.fastq.gz \
    	       -o . \
    	       --taxonomic-profile ${metaphlan_profile_tsv} \
    	       --output-basename ${sample_id} \
    	       --threads 8

    }
    
    output {
    	File humann_out_genefamilies = "${sample_id}_genefamilies.tsv"
    	File humann_out_pathabundance = "${sample_id}_pathabundance.tsv"
    	File humann_out_pathcoverage = "${sample_id}_pathcoverage.tsv"
    	File humann_log = "${sample_id}_humann_temp/${sample_id}.log"
    }

	runtime {
		docker: "gcr.io/osullivan-lab/humann:v3.6"
		cpu: 8 
	  	memory: if input_file_size > 3 then "48GB" else "32GB" 
  		disks: "local-disk 250 HDD"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
	}
}
