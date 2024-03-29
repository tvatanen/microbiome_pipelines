workflow workflowBiobakery {

	File Sample_Path_List
	String Fastq1_Extension
	String Fastq2_Extension
	Map[String, String] SamplesPaths = read_map(Sample_Path_List)
	scatter (pair in SamplesPaths){

		String SampleDir = pair.right
		String sample = pair.left

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

		call metaphlan {
			input: 
			sample=sample, 
			r1=qcQualityHuman.fileR1, 
			r2=qcQualityHuman.fileR2, 
			s1=qcQualityHuman.fileS1, 
			s2=qcQualityHuman.fileS2,
			sample=sample		
		}
		
		call humann2 {
			input: 
			sample=sample, 
			r1=qcQualityHuman.fileR1, 
			r2=qcQualityHuman.fileR2, 
			s1=qcQualityHuman.fileS1, 
			s2=qcQualityHuman.fileS2,
			taxProfile=metaphlan.fileProfile,
			sample=sample		
		}

		call regroupHumann2 {
			input:
			geneFamilies=humann2.fileGeneFamilies,
			sample=sample
		}

		call strainphlanMarkers {
			input:
			sampleBAM=metaphlan.fileBam,
			sample=sample
		}
	}

	call kneaddataReadCountTable {
		input:
		logFiles=qcQualityHuman.log_file
	}

	call combineMetaphlan {
		input:
		taxProfiles=metaphlan.fileProfile
	}

	call combineHumann2 {
		input:
		filesGeneFamilies=humann2.fileGeneFamilies,
		filesPathways=humann2.filePathwayAbundance,
		filesEC=regroupHumann2.fileECAbundance,
		filesKO=regroupHumann2.fileKOAbundance,
		filesEggnog=regroupHumann2.fileEggnogAbundance,
		filesGO1000=regroupHumann2.fileGO1000Abundance,
		filesPfam=regroupHumann2.filePfamAbundance

	}

	call strainphlanClades {
		input:
		filesSampleMarkers=strainphlanMarkers.sampleMarker
	}

	Array[String] Clades = read_lines(strainphlanClades.cladeList)
	scatter (clade in Clades){
		
		call strainphlanTree {
			input:
			filesSampleMarkers=strainphlanMarkers.sampleMarker,
			clade=clade
		}
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

	command {
		kneaddata --input ${file1} --input ${file2} -o . \
		-db tools-tvat-us/DATABASES/HG19 --trimmomatic-options "HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:50" -t 4 --log ${sample}.log
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
  		disks: "local-disk 501 HDD"
	}
}

task metaphlan {
	File r1
	File r2
	File s1
	File s2
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

	command {
		ls -l /cromwell_root/tools-tvat-us/DATABASES/METAPHLAN2/db_v20/
    	
		zcat ${r1} ${r2} ${s1} ${s2} | metaphlan2.py --bowtie2db /cromwell_root/tools-tvat-us/DATABASES/METAPHLAN2/db_v20 --index v20_m200 --mpa_pkl ${ref1} --input_type multifastq -t rel_ab  --bt2_ps "very-sensitive" --tmp_dir /tmp --ignore_viruses -s ${sample}.sam --bowtie2out ${sample}.bowtie2.out > ${sample}.relative_abundance.txt

        # Copy bam file
        samtools view -bS ${sample}.sam -o ${sample}.bam
    }
    
    output {
    	File fileProfile = "${sample}.relative_abundance.txt"
    	File fileBam = "${sample}.bam"
    	File fileBowtie = "${sample}.bowtie2.out"
    }
	
	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "8GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 50 HDD"
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
    	tar -xf ${refChocophlan}

    	#prepare fastqs
    	zcat ${r1} ${r2} ${s1} ${s2} > ${sample}.fq
    	#prepare output folder
    	mkdir humann2_run
		humann2 -i ${sample}.fq -o humann2_run --taxonomic-profile ${taxProfile} --threads 8 --translated-alignment diamond --search-mode uniref90 --gap-fill on --nucleotide-database chocophlan --protein-database /cromwell_root/tools-tvat-us/DATABASES/HUMANN2/UNIREF --bowtie2 /appdownload/bowtie2-2.3.4.1-linux-x86_64 
    }
    
    output {
    	File fileGeneFamilies = "humann2_run/${sample}_genefamilies.tsv"
    	File filePathwayAbundance = "humann2_run/${sample}_pathabundance.tsv"
    	File filePathwayCoverage = "humann2_run/${sample}_pathcoverage.tsv"
    	File fileLog = "humann2_run/${sample}_humann2_temp/${sample}.log"
    }
	
	runtime {
		docker:"gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 8 
	  	memory: if input_file_size > 3 then "64GB" else "24GB" 
  		disks: "local-disk 250 HDD"
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

task combineMetaphlan {
	Array[File] taxProfiles

	command {

	    cat ${write_lines(taxProfiles)} > taxProfiles_2_download.txt
        mkdir dir_taxProfiles
        cat taxProfiles_2_download.txt | gsutil -m cp -I dir_taxProfiles/
		
        /appdownload/metaphlan2/utils/merge_metaphlan_tables.py dir_taxProfiles/* > metaphlan_merged_results.tsv
	}

	output {
		File metaphlanMerged = "metaphlan_merged_results.tsv"
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
	Array[File] filesGeneFamilies
	Array[File] filesPathways
	Array[File] filesEC
	Array[File] filesKO
	Array[File] filesEggnog
	Array[File] filesGO1000
	Array[File] filesPfam

	command {
        cat ${write_lines(filesGeneFamilies)} > filesGeneFamilies_2_download.txt
        mkdir dir_filesGeneFamilies
        cat filesGeneFamilies_2_download.txt | gsutil -m cp -I dir_filesGeneFamilies/
        tar -zcf genefamilies.tar.gz dir_filesGeneFamilies
        rm -r dir_filesGeneFamilies/

        cat ${write_lines(filesPathways)} > filesPathways_2_download.txt
        mkdir dir_filesPathways
        cat filesPathways_2_download.txt | gsutil -m cp -I dir_filesPathways/
        humann2_join_tables -i dir_filesPathways -o pathways_merged_results.tsv
		rm -r dir_filesPathways/
        
        cat ${write_lines(filesEC)} > filesEC_2_download.txt
        mkdir dir_filesEC
        cat filesEC_2_download.txt | gsutil -m cp -I dir_filesEC/
        humann2_join_tables -i dir_filesEC -o ec_merged_results.tsv
		rm -r dir_filesEC/
        
        cat ${write_lines(filesKO)} > filesKO_2_download.txt
        mkdir dir_filesKO
        cat filesKO_2_download.txt | gsutil -m cp -I dir_filesKO/
        humann2_join_tables -i dir_filesKO -o ko_merged_results.tsv
		rm -r dir_filesKO/
        
        cat ${write_lines(filesEggnog)} > filesEggnog_2_download.txt
        mkdir dir_filesEggnog
        cat filesEggnog_2_download.txt | gsutil -m cp -I dir_filesEggnog/
        humann2_join_tables -i dir_filesEggnog -o eggnog_merged_results.tsv
		rm -r dir_filesEggnog/
        
        cat ${write_lines(filesGO1000)} > filesGO1000_2_download.txt
        mkdir dir_filesGO1000
        cat filesGO1000_2_download.txt | gsutil -m cp -I dir_filesGO1000/
        humann2_join_tables -i dir_filesGO1000 -o go1000_merged_results.tsv
		rm -r dir_filesGO1000/
        
        cat ${write_lines(filesPfam)} > filesPfam_2_download.txt
        mkdir dir_filesPfam
        cat filesPfam_2_download.txt | gsutil -m cp -I dir_filesPfam/
        humann2_join_tables -i dir_filesPfam -o pfam_merged_results.tsv
		rm -r dir_filesPfam/

	}

	output {
		File genefamiliesMerged = "genefamilies.tar.gz"
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

task strainphlanMarkers {
	
	File sampleBAM
	String sample
	
	command {
		samtools-0.1.19 view -h ${sampleBAM} -o ${sample}.sam
		/appdownload/metaphlan2/strainphlan_src/sample2markers.py --min_read_depth 5 --ifn_samples ${sample}.sam --input_type sam --output_dir . --nprocs 2 --verbose --samtools_exe /app/samtools-0.1.19
	}

	output {
		
		File sampleMarker = "${sample}.markers"
	}

	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 2
  		memory: "10GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 100 HDD"
	}

}

task strainphlanClades {

	Array [File] filesSampleMarkers
	File refPKL

	command <<<
		/appdownload/metaphlan2/strainphlan.py --mpa_pkl ${refPKL} --ifn_samples ${sep=' ' filesSampleMarkers} --output_dir . --nprocs_main 1 --print_clades_only > clades.txt
		cat clades.txt | tr -d "'" | sed 's/,.*//' | sed 's/(//' > clades_edited.txt
	>>>

	output {
		File cladeList = "clades_edited.txt"
	}

	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "10GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 100 HDD"
	}

}

task strainphlanTree {
	String clade
	File refPKL
	File refAllMarkers
	File refStrainPhlanDB
	Array [File] filesSampleMarkers

	command <<<
		# get the strainphlan database
    	tar -zxf ${refStrainPhlanDB} --strip-components=3
		
		/appdownload/metaphlan2/strainphlan_src/extract_markers.py --mpa_pkl ${refPKL} --ifn_markers ${refAllMarkers} --clade ${clade} --ofn_markers ${clade}.markers.fasta

		mkdir ${clade}

    	if [ -d strainphlan_ref/${clade} ]; then
			/appdownload/metaphlan2/strainphlan.py --mpa_pkl ${refPKL} --ifn_samples ${sep=' ' filesSampleMarkers} --ifn_markers ${clade}.markers.fasta --ifn_ref_genomes strainphlan_ref/${clade}/* --output_dir ${clade} --nprocs_main 16 --clades ${clade} --relaxed_parameters3
		else
			/appdownload/metaphlan2/strainphlan.py --mpa_pkl ${refPKL} --ifn_samples ${sep=' ' filesSampleMarkers} --ifn_markers ${clade}.markers.fasta --output_dir ${clade} --nprocs_main 16 --clades ${clade} --relaxed_parameters3
		fi

		tar -czf ${clade}.tar.gz ${clade}
		
    >>>
	
	output {
		File strainphlanTar = "${clade}.tar.gz"
	}

	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 16
  		memory: "104GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 100 HDD"
	}
}
