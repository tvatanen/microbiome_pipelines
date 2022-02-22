workflow kneaddata_summary_table {

	call run_kneaddata_read_count_table

}

task run_kneaddata_read_count_table {
    Array[File] logFiles

    command {

        mkdir dir_logfiles
        cat ${write_lines(logFiles)} > logfiles.txt
        while read log_file; do
            cp $log_file dir_logfiles
        done <logfiles.txt

        kneaddata_read_count_table --input dir_logfiles/ --output kneaddata_read_count_table.tsv

    }

    output {
        File kneaddata_read_count_table = "kneaddata_read_count_table.tsv"
    }

    runtime {
        docker: "gcr.io/osullivan-lab/kneaddata:v0.10.0"
        cpu: 1
        memory: "4GB"
        preemptible: 2
        disks: "local-disk 50 HDD"
        maxRetries: 2
    }
}