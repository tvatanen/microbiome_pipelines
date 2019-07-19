workflow workflowBiobakery {

    File Sample_Path_List
    String Fastq_Extension
    Map[String, String] SamplesPaths = read_map(Sample_Path_List)
    scatter (pair in SamplesPaths){

		String SampleDir = pair.right
		String sample = pair.left

        File input_file = SampleDir + sample + Fastq_Extension
		call qcAdapters {
			input: 
			sample=sample, 
			file=input_file, 
		}
		
		call qcQualityHuman {
			input: 
			sample=sample, 
			file=qcAdapters.output_file, 
		}

		call metaphlan {
			input: 
			sample=sample, 
			file=qcQualityHuman.output_file, 
			sample=sample		
		}
		
		call humann2 {
			input: 
			sample=sample, 
			file=qcQualityHuman.output_file, 
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
		filesPathways=humann2.filePathwayAbundamce,
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
    File file
    String sample

    command {

        # move file into name that fits with the sample naming strategy
        mv ${file} ${sample}.fq.gz

        trim_galore --phred33 --quality 0 --stringency 5 --length 10 \
        ${sample}.fq.gz

        mv ${sample}_trimmed.fq.gz ${sample}.adapterTrimmed.fq.gz
    }
    
    output {
        File output_file = "${sample}.adapterTrimmed.fq.gz"
    }
	
	runtime {
		docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "1GB"
  		preemptible: 2
  		disks: "local-disk 40 HDD"
	}
}

task qcQualityHuman {
	File file
	String sample
	File ref1
	File ref2
	File ref3
	File ref4
	File ref5
	File ref6

    command {
        kneaddata --input ${file} -o . \
        -db $(dirname ${ref1}) --trimmomatic-options "HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:50" -t 4 --log ${sample}_kneaddata.log
        rm *trimmed*
        rm *bowtie2*
        gzip ${sample}.adapterTrimmed_kneaddata.fastq
    }

    output {
        File output_file = "${sample}.adapterTrimmed_kneaddata.fastq.gz"
        File log_file = "${sample}_kneaddata.log"
    }
	
	runtime {
		docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 4
  		memory: "24GB"
  		preemptible: 2
  		disks: "local-disk 501 HDD"
	}
}

task metaphlan {
	File file
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
    	ls -l /cromwell_root/tools-tvat/DATABASES/METAPHLAN2/db_v20/
    	
    	zcat ${file} | metaphlan2.py --bowtie2db /cromwell_root/tools-tvat/DATABASES/METAPHLAN2/db_v20 --index v20_m200 --mpa_pkl ${ref1} --input_type multifastq -t rel_ab  --bt2_ps "very-sensitive" --tmp_dir /tmp --ignore_viruses -s ${sample}.sam --bowtie2out ${sample}.bowtie2.out > ${sample}.relative_abundance.txt

        # Copy bam file
        samtools view -bS ${sample}.sam -o ${sample}.bam
    }
    
    output {
    	File fileProfile = "${sample}.relative_abundance.txt"
    	File fileBam = "${sample}.bam"
    	File fileBowtie = "${sample}.bowtie2.out"
    }
	
	runtime {
		docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "8GB"
  		preemptible: 2
  		disks: "local-disk 50 HDD"
	}

}

task humann2 {
	File file
	File taxProfile
	File refChocophlan
	File refUniref90
	String sample
	
	command {    	
    	# get the chocophlan database
    	tar -zxvf ${refChocophlan}

    	#prepare fastqs
    	zcat ${file} > ${sample}.fq
    	#prepare output folder
    	mkdir humann2_run
    	humann2 -i ${sample}.fq -o humann2_run --taxonomic-profile ${taxProfile} --threads 8 --translated-alignment diamond --search-mode uniref90 --gap-fill on --nucleotide-database chocophlan --protein-database /cromwell_root/tools-tvat/DATABASES/HUMANN2/UNIREF --bowtie2 /appdownload/bowtie2-2.3.4.1-linux-x86_64 
    }
    
    output {
    	File fileGeneFamilies = "humann2_run/${sample}_genefamilies.tsv"
    	File filePathwayAbundamce = "humann2_run/${sample}_pathabundance.tsv"
    	File filePathwayCoverage = "humann2_run/${sample}_pathcoverage.tsv"
    	File fileLog = "humann2_run/${sample}_humann2_temp/${sample}.log"
    }
	
	runtime {
		docker:"asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 8
  		memory: "24GB"
  		disks: "local-disk 120 HDD"
  		preemptible: 2
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
        docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
        cpu: 1
        memory: "4GB"
        preemptible: 2
        disks: "local-disk 50 HDD"
    }

}

task kneaddataReadCountTable {
    Array[File] logFiles

    command {
        mkdir logfiles

        for logfile in $(cat ${write_lines(logFiles)})
        do
            cp $logfile logfiles
        done

        kneaddata_read_count_table --input logfiles/ --output kneaddata_read_count_table.tsv

    }

    output {
        File kneaddataReadCountTable = "kneaddata_read_count_table.tsv"
    }

    runtime {
        docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
        cpu: 1
        memory: "4GB"
        preemptible: 2
        disks: "local-disk 50 HDD"
    }
}

task combineMetaphlan {
    Array[File] taxProfiles

    command {
        /appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " taxProfiles} > metaphlan_merged_results.tsv
    }

    output {
        File metaphlanMerged = "metaphlan_merged_results.tsv"
    }

    runtime {
        docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
        cpu: 1
        memory: "4GB"
        preemptible: 2
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
        mkdir genefamilies
        cp -t genefamilies ${sep=" " filesGeneFamilies}
        humann2_join_tables -i genefamilies -o genefamilies_merged_results.tsv

        mkdir pathways
        cp -t pathways ${sep=" " filesPathways}
        humann2_join_tables -i pathways -o pathways_merged_results.tsv

        mkdir ec
        cp -t ec ${sep=" " filesEC}
        humann2_join_tables -i pathways -o ec_merged_results.tsv

        mkdir ko
        cp -t ko ${sep=" " filesKO}
        humann2_join_tables -i ko -o ko_merged_results.tsv

        mkdir eggnog
        cp -t eggnog ${sep=" " filesEggnog}
        humann2_join_tables -i eggnog -o eggnog_merged_results.tsv

        mkdir go
        cp -t go ${sep=" " filesGO1000}
        humann2_join_tables -i go -o go1000_merged_results.tsv

        mkdir pfam
        cp -t pfam ${sep=" " filesPfam}
        humann2_join_tables -i pfam -o pfam_merged_results.tsv

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
        docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
        cpu: 1
        memory: "100GB"
        preemptible: 2
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
        docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
        cpu: 2
        memory: "10GB"
        preemptible: 2
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
        docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
        cpu: 1
        memory: "10GB"
        preemptible: 2
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
        tar -zxvf ${refStrainPhlanDB} --strip-components=3
        
        /appdownload/metaphlan2/strainphlan_src/extract_markers.py --mpa_pkl ${refPKL} --ifn_markers ${refAllMarkers} --clade ${clade} --ofn_markers ${clade}.markers.fasta

        mkdir ${clade}

        if [ -d strainphlan_ref/${clade} ]; then
            /appdownload/metaphlan2/strainphlan.py --mpa_pkl ${refPKL} --ifn_samples ${sep=' ' filesSampleMarkers} --ifn_markers ${clade}.markers.fasta --ifn_ref_genomes strainphlan_ref/${clade}/* --output_dir ${clade} --nprocs_main 16 --clades ${clade}
        else
            /appdownload/metaphlan2/strainphlan.py --mpa_pkl ${refPKL} --ifn_samples ${sep=' ' filesSampleMarkers} --ifn_markers ${clade}.markers.fasta --output_dir ${clade} --nprocs_main 16 --clades ${clade}
        fi

        tar -czf ${clade}.tar.gz ${clade}
        
    >>>
    
    output {
        File strainphlanTar = "${clade}.tar.gz"
    }

    runtime {
        docker: "asia.gcr.io/osullivan-lab/metagenomicstools:03072019"
        cpu: 16
        memory: "104GB"
        preemptible: 2
        disks: "local-disk 100 HDD"
    }
}
