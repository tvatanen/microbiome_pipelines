workflow strainphlan_markers {
	
	call sample2markers

}

task sample2markers {
	File sampleBAM
	String sample_id

	command {

		samtools view -h ${sampleBAM} -o ${sample_id}.sam
		bzip2 ${sample_id}.sam
		sample2markers.py -i ${sample_id}.sam.bz2 -o . --nproc 1

	}

	output {
		
		File strainphlan_markers = "${sample_id}.pkl"
	}

	runtime {
		docker: "gcr.io/osullivan-lab/metaphlan:v3.0.12"
		cpu: 1
  		memory: "8GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 50 HDD"
	}
}
