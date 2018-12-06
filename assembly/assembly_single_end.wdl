workflow workflowAssembly {

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
            input_file=input_file, 
        }
        call qcQualityHuman {
            input: 
            sample=sample, 
            input_file=qcAdapters.output_file, 
        }
        call assemble {
            input: 
            sample=sample, 
            input_file=qcQualityHuman.output_file, 
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

    scatter (pair in SamplesPaths){
        String SampleDir_map = pair.right
        String sample_map = pair.left

        File input_file_map = SampleDir_map + sample_map + Fastq_Extension

        call map { 
            input: 
            input_file=input_file_map,
            sample=sample_map,
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
    File input_file
    String sample

    command {

        # move file into name that fits with the sample naming strategy
        mv ${input_file} ${sample}.fq.gz

        trim_galore --phred33 --quality 0 --stringency 5 --length 10 \
        ${sample}.fq.gz

        mv ${sample}_trimmed.fq.gz ${sample}.adapterTrimmed.fq.gz
    }
    
    output {
        File output_file = "${sample}.adapterTrimmed.fq.gz"
    }
    
    runtime {
        docker: "tvatanen/metagenomictools"
    }
}

task qcQualityHuman {
    File input_file
    String sample
    File ref1
    File ref2
    File ref3
    File ref4
    File ref5
    File ref6

    command {
        kneaddata --input ${input_file} -o . \
        -db $(dirname ${ref1}) --trimmomatic-options "HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:50" -t 4
        rm *trimmed*
        rm *bowtie2*
        gzip ${sample}.adapterTrimmed_kneaddata.fastq
    }

    output {
        File output_file = "${sample}.adapterTrimmed_kneaddata.fastq.gz"
    }
        
    runtime {
        docker: "tvatanen/metagenomictools:java1.8"
    }
}

task assemble {
    File input_file
    String sample

    command <<<
        rm -f assemble
        megahit -r ${input_file} -t 4 -m 15000000000 -o assemble
        cat assemble/final.contigs.fa | \
        awk -v var="${sample}" '
            {if($0 ~ /^>/) {contigName=substr($0, 2,length($0))} 
            else {seq=$0; if(length($0) >= 500) {print ">"var"_"contigName"\n"seq}} }' > assemble/${sample}.min500.contigs.fa
            >>>
    output {
        File fileContigs = "assemble/${sample}.min500.contigs.fa"
    }
    runtime {
        docker: "tvatanen/metagenomictools"
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
        docker: "tvatanen/metagenomictools"
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
        docker: "tvatanen/metagenomictools"
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
        docker: "tvatanen/metagenomictools"
    }
}

task map {
    File input_file
    String sample
    File nrFa
    File nrFai
    File ref1
    File ref2
    File ref3
    File ref4
    File ref5

    command {
        
        bwa mem -t 8 -M ${nrFa} ${input_file} | \
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
        docker: "tvatanen/metagenomictools"
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
        docker: "tvatanen/metagenomictools"
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
        docker: "tvatanen/r-docker:v3.4.4"
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
        docker: "tvatanen/r-docker:v3.4.4"
    }
}





