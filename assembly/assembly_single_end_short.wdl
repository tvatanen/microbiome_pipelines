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