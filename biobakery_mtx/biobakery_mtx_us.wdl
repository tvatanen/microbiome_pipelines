workflow workflowMTXBiobakery {

	File Sample_Path_List
	String Fastq1_Extension
	String Fastq2_Extension
	Array[Array[String]] SamplesPaths = read_tsv(Sample_Path_List)
	
	scatter (pair in SamplesPaths){
		String SampleDir = pair[1]
		String sample = pair[0]
		File MetaphlanProfile = pair[2]

		File F1 = SampleDir + sample + Fastq1_Extension
		File F2 = SampleDir + sample + Fastq2_Extension

		call qcAdapters {
			input: 
			sample=sample, 
			file1=F1, 
			file2=F2
		}
		
		call qcQualityHuman {
			input: 
			sample=sample, 
			file1=qcAdapters.fileR1, 
			file2=qcAdapters.fileR2
		}
		
		call humann2 {
			input: 
			sample=sample, 
			r1=qcQualityHuman.fileR1, 
			r2=qcQualityHuman.fileR2, 
			s1=qcQualityHuman.fileS1, 
			s2=qcQualityHuman.fileS2,
			taxProfile=MetaphlanProfile,
			sample=sample		
		}

		call regroupHumann2 {
			input:
			geneFamilies=humann2.fileGeneFamilies,
			sample=sample
		}

	}

	call kneaddataReadCountTable {
		input:
		logFiles=qcQualityHuman.log_file
	}

	call combineHumann2 {
		input:
		filesGeneFamilies=humann2.fileGeneFamilies,
		filesPathways=humann2.filePathwayAbundamce,
		filesEC=regroupHumann2.fileECAbundance,
		filesKO=regroupHumann2.fileKOAbundance,
		filesEggnog=regroupHumann2.fileEggnogAbundance,
		filesGO1000=regroupHumann2.fileGO1000Abundance,
		filesPfam=regroupHumann2.filePfamAbundance

	}

}


task qcAdapters {
	File file1
	File file2
	String sample

	command {

		# move file into name that fits with the sample naming strategy
		mv ${file1} ${sample}.1.fq.gz
		mv ${file2} ${sample}.2.fq.gz

		trim_galore --paired --phred33 --quality 0 --stringency 5 --length 10 \
		${sample}.1.fq.gz ${sample}.2.fq.gz

		mv ${sample}.1_val_1.fq.gz ${sample}.adapterTrimmed.1.fq.gz
		mv ${sample}.2_val_2.fq.gz ${sample}.adapterTrimmed.2.fq.gz
	}
	
	output {
		File fileR1 = "${sample}.adapterTrimmed.1.fq.gz"
		File fileR2 = "${sample}.adapterTrimmed.2.fq.gz"

	}
	
	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "1GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
        disks: "local-disk 40 HDD"
	}
}


task qcQualityHuman {
	File file1
	File file2
	String sample
	File ref1
	File ref2
	File ref3
	File ref4
	File ref5
	File ref6
	File ref7
	File ref8
	File ref9
	File ref10
	File ref11
	File ref12
	File ref13
	File ref14
	File ref15
	File ref16
	File ref17
	File ref18

	command {
		kneaddata --input ${file1} --input ${file2} -o . \
		-db tools-tvat-us/DATABASES/HG19/Homo_sapiens \
		-db tools-tvat-us/DATABASES/HOST_TRANSCRIPT/rna-bowtie2 \
		-db tools-tvat-us/DATABASES/SILVA_Bowtie/rna.fna \
		--trimmomatic-options "HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:50" -t 4 \
		--log ${sample}.log
		rm *trimmed*
		rm *bowtie2*
		
		gzip ${sample}.adapterTrimmed.1_kneaddata_paired_1.fastq
		gzip ${sample}.adapterTrimmed.1_kneaddata_paired_2.fastq
		gzip ${sample}.adapterTrimmed.1_kneaddata_unmatched_1.fastq
		gzip ${sample}.adapterTrimmed.1_kneaddata_unmatched_2.fastq
	}
	
	output {
		File fileR1 = "${sample}.adapterTrimmed.1_kneaddata_paired_1.fastq.gz"
		File fileR2 = "${sample}.adapterTrimmed.1_kneaddata_paired_2.fastq.gz"
		File fileS1 = "${sample}.adapterTrimmed.1_kneaddata_unmatched_1.fastq.gz"
		File fileS2 = "${sample}.adapterTrimmed.1_kneaddata_unmatched_2.fastq.gz"
		File log_file = "${sample}.log"
	}
	
	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 4
  		memory: "24GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 1000 HDD"
	}
}

task humann2 {
	File r1
	File r2
	File s1
	File s2
	File taxProfile
	File refChocophlan
	File refUniref90
	String sample
	Float input_file_size = size(r1, "G") + size(r2, "G") + size(s1, "G") + size(s2, "G")
	
	command {    	
    	# get the chocophlan database
    	tar -zxf ${refChocophlan}

    	#prepare fastqs
    	zcat ${r1} ${r2} ${s1} ${s2} > ${sample}.fq
    	#prepare output folder
    	mkdir humann2_run
    	humann2 -i ${sample}.fq -o humann2_run --taxonomic-profile ${taxProfile} --threads 8 --translated-alignment diamond --search-mode uniref90 --gap-fill on --nucleotide-database chocophlan --protein-database /cromwell_root/tools-tvat-us/DATABASES/HUMANN2/UNIREF --bowtie2 /appdownload/bowtie2-2.3.4.1-linux-x86_64 
    }
    
    output {
    	File fileGeneFamilies = "humann2_run/${sample}_genefamilies.tsv"
    	File filePathwayAbundamce = "humann2_run/${sample}_pathabundance.tsv"
    	File filePathwayCoverage = "humann2_run/${sample}_pathcoverage.tsv"
    	File fileLog = "humann2_run/${sample}_humann2_temp/${sample}.log"
    }

	runtime {
		docker:"gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 8 
	  	memory: if input_file_size > 3 then "64GB" else "24GB" 
  		disks: "local-disk 500 HDD"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
	}
	
}

task regroupHumann2 {
	File geneFamilies
	File refEC
	File refKO
	File refPfam
	File refEggnog
	File refGO1000
	String sample

	command {
		humann2_regroup_table -i ${geneFamilies} -c ${refEC} -o ${sample}_ec.tsv
		humann2_regroup_table -i ${geneFamilies} -c ${refKO} -o ${sample}_ko.tsv
		humann2_regroup_table -i ${geneFamilies} -c ${refPfam} -o ${sample}_pfam.tsv
		humann2_regroup_table -i ${geneFamilies} -c ${refEggnog} -o ${sample}_eggnog.tsv
		humann2_regroup_table -i ${geneFamilies} -c ${refGO1000} -o ${sample}_go1000.tsv

	}

	output {
		File fileECAbundance = "${sample}_ec.tsv"
		File fileKOAbundance = "${sample}_ko.tsv"
		File filePfamAbundance = "${sample}_pfam.tsv"
		File fileEggnogAbundance = "${sample}_eggnog.tsv"
		File fileGO1000Abundance = "${sample}_go1000.tsv"

	}	

	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "4GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 50 HDD"
	}
}

task kneaddataReadCountTable {
	Array[File] logFiles

	command {

		cat ${write_lines(logFiles)} > logFiles_2_download.txt
        mkdir dir_logfiles
        cat logFiles_2_download.txt | gsutil -m cp -I dir_logfiles/
        
		kneaddata_read_count_table --input dir_logfiles/ --output kneaddata_read_count_table.tsv

	}

	output {
		File kneaddataReadCountTable = "kneaddata_read_count_table.tsv"
	}

	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "4GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 50 HDD"
	}
}


task combineHumann2 {
	Array[String] filesGeneFamilies
	Array[String] filesPathways
	Array[String] filesEC
	Array[String] filesKO
	Array[String] filesEggnog
	Array[String] filesGO1000
	Array[String] filesPfam

	command {

		cat ${write_lines(filesGeneFamilies)} > filesGeneFamilies_2_download.txt
		mkdir dir_filesGeneFamilies
		cat filesGeneFamilies_2_download.txt | gsutil -m cp -I dir_filesGeneFamilies
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py dir_filesGeneFamilies/* > genefamilies_merged_results.tsv

		cat ${write_lines(filesPathways)} > filesPathways_2_download.txt
		mkdir dir_filesPathways
		cat filesPathways_2_download.txt | gsutil -m cp -I dir_filesPathways
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py dir_filesPathways/* > pathways_merged_results.tsv
		
		cat ${write_lines(filesEC)} > filesEC_2_download.txt
		mkdir dir_filesEC
		cat filesEC_2_download.txt | gsutil -m cp -I dir_filesEC
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py dir_filesEC/* > ec_merged_results.tsv

		cat ${write_lines(filesKO)} > filesKO_2_download.txt
		mkdir dir_filesKO
		cat filesKO_2_download.txt | gsutil -m cp -I dir_filesKO
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py dir_filesKO/* > ko_merged_results.tsv

		cat ${write_lines(filesEggnog)} > filesEggnog_2_download.txt
		mkdir dir_filesEggnog
		cat filesEggnog_2_download.txt | gsutil -m cp -I dir_filesEggnog
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py dir_filesEggnog/* > eggnog_merged_results.tsv

		cat ${write_lines(filesGO1000)} > filesGO1000_2_download.txt
		mkdir dir_filesGO1000
		cat filesGO1000_2_download.txt | gsutil -m cp -I dir_filesGO1000
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py dir_filesGO1000/* > go1000_merged_results.tsv
		
		cat ${write_lines(filesPfam)} > filesPfam_2_download.txt
		mkdir dir_filesPfam
		cat filesPfam_2_download.txt | gsutil -m cp -I dir_filesPfam
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py dir_filesPfam/* > pfam_merged_results.tsv

	}

	output {
		File genefamiliesMerged = "genefamilies_merged_results.tsv"
		File pathwaysMerged = "pathways_merged_results.tsv"
		File ecMerged = "ec_merged_results.tsv"
		File koMerged = "ko_merged_results.tsv"
		File eggnogMerged = "eggnog_merged_results.tsv"
		File go1000Merged = "go1000_merged_results.tsv"
		File pfamMerged = "pfam_merged_results.tsv"

	}

	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "100GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 500 HDD"
	}

}