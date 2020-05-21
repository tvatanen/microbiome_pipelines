workflow assembly_binning {

    File Sample_Path_List
    String Fastq1_Extension
    String Fastq2_Extension
    Array[Array[String]] SamplesPaths = read_tsv(Sample_Path_List)
    scatter (pair in SamplesPaths){

        String SampleDir = pair[1]
        String sample = pair[0]
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
        call map { 
            input: 
            fileR1=qcQualityHuman.fileR1,
            fileR2=qcQualityHuman.fileR2,
            sample=sample,
            contigs=assemble.fileContigs
        }

        call metabat2 {
            input:
            contigs=assemble.fileContigs,
            sample=sample,
            bam=map.fileBAM
        }

        call checkm {
            input:
            bins=metabat2.bins,
            sample=sample
        }

        call gtdbtk {
            input:
            bins=metabat2.bins,
            sample=sample
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
        disks: "local-disk 40 SSD"
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
        disks: "local-disk 100 SSD"
    }
}

task map {
    File fileR1
    File fileR2
    String sample
    File contigs

    command {
        
        bwa index ${contigs}
        samtools faidx ${contigs}

        bwa mem -t 8 -M ${contigs} ${fileR1} ${fileR2} | \
        samtools view - -h -Su -F 2308 -q 0 | \
        samtools sort -@ 8 -m 2G -O bam -o ${sample}.sort.bam 

    }
    
    output {

        File fileBAM = "${sample}.sort.bam"

    }
    
    runtime {
        docker: "gcr.io/osullivan-lab/metagenomicstools:03072019"
        cpu: 8
        memory: "24GB"
        preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
        bootDiskSizeGb: 50
        disks: "local-disk 200 SSD"
    }
}

task metabat2 {
    File bam
    String sample
    File contigs

    command {
        
        runMetaBat.sh ${contigs} ${bam}
        mv ${sample}.min500.contigs.fa.depth.txt ${sample}.depths.txt
        mv ${sample}.min500.contigs.fa*/ ${sample}_bins
        tar -czf ${sample}.bins.tar.gz ${sample}_bins

    }
    
    output {

        File bins = "${sample}.bins.tar.gz"
        File depths = "${sample}.depths.txt"

    }
    
    runtime {
        docker: "gcr.io/osullivan-lab/metabat2:v2.15-3"
        cpu: 8
        memory: "12GB"
        preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
        bootDiskSizeGb: 50
        disks: "local-disk 100 SSD"
    }
}

task checkm {
    File bins
    String sample

    command {

        checkm data setRoot /checkm-data
        tar -xf ${bins}
        checkm lineage_wf -f ${sample}_checkm.txt -t 4 -x fa ${sample}_bins/ ${sample}_checkm
        tar -czf ${sample}_checkm.tar.gz ${sample}_checkm

    }

    output {

        File checkm_summary = "${sample}_checkm.txt"
        File checkm_out = "${sample}_checkm.tar.gz"

    }

    runtime {
        docker: "gcr.io/osullivan-lab/checkm:v1.1.2"
        cpu: 4
        memory: "100GB"
        preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
        disks: "local-disk 100 HDD"
    }
}

task gtdbtk {
    File bins
    File gtdb_reference
    String gtdb_release
    String sample

    command {

        tar -xf ${gtdb_reference} -C /gtdbtk-data/
        export GTDBTK_DATA_PATH=/gtdbtk-data/${gtdb_release}/
        tar -xf ${bins}
        gtdbtk classify_wf --genome_dir ${sample}_bins/ -x fa --cpus 4 --out_dir ${sample}_gtdb
        tar -czf ${sample}_gtdb.tar.gz ${sample}_gtdb
        cp ${sample}_gtdb/classify/gtdbtk.bac120.summary.tsv ${sample}.bac120.summary.tsv

    }

    output {

        File gtdbtk_out = "${sample}_gtdb.tar.gz"
        File gtdbtk_summary = "${sample}.bac120.summary.tsv"

    }

    runtime {
        docker: "gcr.io/osullivan-lab/gtdbtk:v1.0.2"
        cpu: 4
        memory: "128GB"
        preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
        bootDiskSizeGb: 50
        disks: "local-disk 100 SSD"
    }
}

