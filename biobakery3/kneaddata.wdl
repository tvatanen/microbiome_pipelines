workflow qc_kneaddata {

    call kneaddata 

}

task kneaddata {
    File file_r1
    File file_r2
    String sample_id

    command {
        kneaddata --input ${file_r1} \
                  --input ${file_r2} \
                  -o . \
                  --output-prefix ${sample_id} \
                  -db /hg37 \
                  --trimmomatic-options "HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:50" \
                  -t 4 \
                  --reorder \
                  --remove-intermediate-output \
                  --trimmomatic /opt/conda/share/trimmomatic \
                  --log ${sample_id}.log
        
        rm *bowtie2*
        
        gzip ${sample_id}_paired_1.fastq
        gzip ${sample_id}_paired_2.fastq
        gzip ${sample_id}_unmatched_1.fastq
        gzip ${sample_id}_unmatched_2.fastq
    }
    
    output {
        File output_paired1 = "${sample_id}_paired_1.fastq.gz"
        File output_paired2 = "${sample_id}_paired_2.fastq.gz"
        File output_unmatched1 = "${sample_id}_unmatched_1.fastq.gz"
        File output_unmatched2 = "${sample_id}_unmatched_2.fastq.gz"
        File log_file = "${sample_id}.log"
    }
    
    runtime {
        docker: "gcr.io/osullivan-lab/kneaddata:v0.10.0"
        cpu: 4
        memory: "24GB"
        preemptible: 2
        maxRetries: 3
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
        disks: "local-disk 500 HDD"
    }
}
