workflow workflowAssembly {

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
		call assemble {
			input: 
			sample=sample, 
			r1=qcQualityHuman.fileR1, 
			r2=qcQualityHuman.fileR2, 
			s1=qcQualityHuman.fileS1, 
			s2=qcQualityHuman.fileS2
		}
		call predictgenes {
			input: 
			sample=sample, 
			fileContigs=assemble.fileContigs
		}
	}

	call mergeGenePredictions { 
		input: 
		genepredictions=predictgenes.fileFNA 
	}

	call cdhit { 
		input: 
		combinedgenepredictions=mergeGenePredictions.out
	}

	Array[Pair[File, File]] fileR1R2 = zip(qcQualityHuman.fileR1, qcQualityHuman.fileR2) 
	
	scatter (pair in fileR1R2){

			String currentSample = basename(pair.left, ".adapterTrimmed.1_kneaddata_paired_1.fastq.gz")

			call map { 
			input: 
			fileR1=pair.left,
			fileR2=pair.right,
			sample=currentSample,
			nrFa=cdhit.nrFa,
			nrFai=cdhit.nrFai,
			ref1=cdhit.nrRef1,
			ref2=cdhit.nrRef2,
			ref3=cdhit.nrRef3,
			ref4=cdhit.nrRef4,
			ref5=cdhit.nrRef5
		}
	}

	File fileGeneNames = map.fileGenes[1]

	call msp { 
		input: 
		geneCountArray=map.fileCount,
		geneAbundanceArray=map.fileAbundance,
		geneNamesArray=fileGeneNames
	}

	scatter (abundanceFile in map.fileAbundance){
		String sampleName = basename(abundanceFile, ".abundance.txt")

		call msp_abundance {
			input:
			fileAbundance=abundanceFile,
			geneNamesArray=fileGeneNames,
			fileMSPs=msp.fileMSPs,
			sample=sampleName
		}
	}

	call msp_matrix {
		input:
		MSPCoreAbundance=msp_abundance.fileMSPCoreAbundance,
		MSPOtherAbundance=msp_abundance.fileMSPOtherAbundance,
		MSPCore=msp_abundance.fileMSPCore[1],
		MSPOther=msp_abundance.fileMSPOther[1]
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

task assemble {
	File r1
	File r2
	File s1
	File s2
	String sample

	command <<<
		rm -f assemble
		megahit -1 ${r1} -2 ${r2} -r ${s1},${s2} -t 4 -m 15000000000 -o assemble
		cat assemble/final.contigs.fa | \
		awk -v var="${sample}" '
			{if($0 ~ /^>/) {contigName=substr($0, 2,length($0))} 
			else {seq=$0; if(length($0) >= 500) {print ">"var"_"contigName"\n"seq}} }' > assemble/${sample}.min500.contigs.fa
			>>>
	output {
		File fileContigs = "assemble/${sample}.min500.contigs.fa"
	}
	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 4
  		memory: "15GB"
  		preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 100 HDD"
	}
}

task predictgenes {
	File fileContigs
	String sample

	command <<<

		mkdir prodigal
		if [[ `wc -l ${fileContigs} | awk '{print $1}'` == "0" ]]; then
			touch prodigal/${sample}.gff
			touch prodigal/${sample}.fna
			touch prodigal/${sample}.faa
 		else
 			prodigal -p meta -i ${fileContigs} -f gff \
			-o prodigal/${sample}.gff \
			-d prodigal/${sample}.fna \
			-a prodigal/${sample}.faa \
 			2> prodigal/prodigal.stderr > prodigal/prodigal.out
 		fi
	>>>
	
	output {
		File fileFNA = "prodigal/${sample}.fna"
		File fileFAA = "prodigal/${sample}.faa"
		File fileGFF = "prodigal/${sample}.gff"
	}

	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "7GB"
  		preemptible: 2
  		maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 100 HDD"
	}
}

task mergeGenePredictions {
	Array[File] genepredictions
	
	command <<<
		cat ${sep=' ' genepredictions} > combined_genepredictions.fna
	>>>

  	output { 
  		File out = "combined_genepredictions.fna" 
  	}

  	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 1
  		memory: "7GB"
  		preemptible: 2
  		maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 500 HDD"
	}
}

task cdhit {
	File combinedgenepredictions

	command {
		usearch8.1.1861_i86linux64 -sortbylength ${combinedgenepredictions} -fastaout combinedgenepredictions.sorted.fna -minseqlength 102
		cd-hit-est -i combinedgenepredictions.sorted.fna -T 32 -aS 0.9 -c 0.95 -M 0 -r 0 -B 0 -d 0 -o nr.fa
		bwa index nr.fa
		samtools faidx nr.fa
	}

  	output { 
  		File nrFa = "nr.fa"
  		File nrFai = "nr.fa.fai"
  		File nrRef1 = "nr.fa.bwt"
  		File nrRef2 = "nr.fa.pac"
  		File nrRef3 = "nr.fa.ann"
  		File nrRef4 = "nr.fa.amb"
  		File nrRef5 = "nr.fa.sa"
  	}

  	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 32
  		memory: "120GB"
  		bootDiskSizeGb: 50
  		maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 501 HDD"
	}
}

task map {
	File fileR1
	File fileR2
	String sample
	File nrFa
	File nrFai
	File ref1
	File ref2
	File ref3
	File ref4
	File ref5

	command {
		
		bwa mem -t 8 -M ${nrFa} ${fileR1} ${fileR2} | \
		samtools view - -h -Su -F 2308 -q 0 | \
 		samtools sort -n -@ 8 -m 2G -O bam -o ${sample}.sort.bam 

		/app/bamfilter_mate_saver.py -i ${sample}.sort.bam -o ${sample}.sort.filtered.ID95.bam -f 0.95

		samtools flagstat ${sample}.sort.filtered.ID95.bam > ${sample}.ID95.flagstat.txt

		/app/bam2counts.py -t 1 -y 1  -s ${sample} -v ${nrFai} -c ${sample}.count.txt -a ${sample}.abundance.txt -l ${sample}.learn.txt -i ${sample}.sort.filtered.ID95.bam

		cut -f1 ${sample}.count.txt > gene.names.txt
		cut -f2 ${sample}.count.txt > tmp.count.txt
		cut -f2 ${sample}.abundance.txt > tmp.abundance.txt
		
		mv tmp.count.txt ${sample}.count.txt
		mv tmp.abundance.txt ${sample}.abundance.txt

	}
	
	output {
		File fileBAM = "${sample}.sort.bam"
		File fileBAMfiltered = "${sample}.sort.filtered.ID95.bam"
		File fileFlagStat = "${sample}.ID95.flagstat.txt"
		File fileLearn = "${sample}.learn.txt"
		File fileCount = "${sample}.count.txt"
		File fileAbundance = "${sample}.abundance.txt"
		File fileGenes = "gene.names.txt"


	}
	
	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 8
  		memory: "24GB"
  		preemptible: 2
  		maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		bootDiskSizeGb: 50
  		disks: "local-disk 200 HDD"
	}
}

task msp {
	Array[File] geneCountArray
	Array[File] geneAbundanceArray
	File geneNamesArray

	command <<<
	
	mkdir msps
	
	paste ${geneNamesArray} ${sep=' ' geneCountArray} > combined_geneCount.txt
	paste ${geneNamesArray} ${sep=' ' geneAbundanceArray} > combined_geneAbundance.txt

	mspminer /appdownload/mspminer-bin/settings.ini
	tar -czf msps.tar.gz msps
	
	>>>
	
	output {
		File fileCountMatrix = "combined_geneCount.txt"
		File fileAbundanceMatrix = "combined_geneAbundance.txt"
		File fileMSPs = "msps.tar.gz"	
	}
	
	runtime {
		docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
		cpu: 32
  		memory: "120GB"
  		bootDiskSizeGb: 50
  		maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 500 HDD"
	}
}

task msp_abundance {
	File fileAbundance
	File geneNamesArray
	String sample
	File fileMSPs
	
	command <<<

		sum_kbp=`tail -n +2 ${fileAbundance} | awk '{ sum+= $1 } END {print sum}'`

		multiplier=`echo | awk -v million_scale=1000000 -v sum_kbp="$sum_kbp" '{printf "%.5g", million_scale/sum_kbp}'`

		head -1 ${fileAbundance} > ${sample}.TPM.txt
		tail -n +2 ${fileAbundance} | awk -v multiplier="$multiplier" '{ printf "%.5g\n", $1*multiplier }' >> ${sample}.TPM.txt
		
		paste ${geneNamesArray} ${sample}.TPM.txt > geneNames.${sample}.TPM.txt

		cp ${fileMSPs} msps.tar.gz
		tar -zxvf msps.tar.gz

		for folder in `ls msps | grep msp`; do
        	echo $folder
        	paste <(yes $folder | head -n `wc -l msps/$folder/modules.tsv | awk '{print $1-1}'`) <(tail -n +2 msps/$folder/modules.tsv) >> combined_msp_gene_table.txt
		done

		Rscript --vanilla --slave /app/get_msp_abundance.R combined_msp_gene_table.txt geneNames.${sample}.TPM.txt ${sample}

		cut -f1 ${sample}.core_msp_30_TPM_median.txt > names_msp_core.txt
		cut -f2 ${sample}.core_msp_30_TPM_median.txt > tmp_core.txt
		mv tmp_core.txt ${sample}.core_msp_30_TPM_median.txt

		cut -f1 ${sample}.other_msp_30_TPM_median.txt > names_msp_other.txt
		cut -f2 ${sample}.other_msp_30_TPM_median.txt > tmp_other.txt
		mv tmp_other.txt ${sample}.other_msp_30_TPM_median.txt
	>>>

	output {
		File fileMSPCoreAbundance = "${sample}.core_msp_30_TPM_median.txt"
		File fileMSPOtherAbundance = "${sample}.other_msp_30_TPM_median.txt"
		File fileMSPCore = "names_msp_core.txt"
		File fileMSPOther = "names_msp_other.txt"
		File fileTPM = "${sample}.TPM.txt"	
	}
	
	runtime {
		docker: "gcr.io/osullivan-lab/r-docker:v3.4.4"
		cpu: 2
  		memory: "10GB"
  		bootDiskSizeGb: 50
  		maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 100 HDD"
	}

}


task msp_matrix {
	Array[File] MSPCoreAbundance
	Array[File] MSPOtherAbundance
	File MSPCore
	File MSPOther

	command <<<
		
	paste ${MSPCore} ${sep=' ' MSPCoreAbundance} > combined_MSPCoreAbundance.txt
	paste ${MSPOther} ${sep=' ' MSPOtherAbundance} > combined_MSPOtherAbundance.txt

	>>>
	
	output {
		File fileMSPCoreAbundance = "combined_MSPCoreAbundance.txt"
		File fileMSPOtherAbundance = "combined_MSPOtherAbundance.txt"
	}
	
	runtime {
		docker: "gcr.io/osullivan-lab/r-docker:v3.4.4"
		cpu: 1
  		memory: "10GB"
  		bootDiskSizeGb: 50
  		maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
  		disks: "local-disk 100 HDD"
	}
}





