workflow annotate_gene_catalogue {
    Array[File] gene_clusters_split
    Int preemptible_tries

    scatter (gene_shard in gene_clusters_split){

        call annotate_eggnog {
            input:
            num_preemptible=preemptible_tries,
            gene_catalogue=gene_shard 
        }
    }
}

task annotate_eggnog {
    File gene_catalogue
    File eggnog_db
    File eggnog_db_diamond
    Int eggnog_memory_gb
    Int eggnog_num_cores
    Int num_preemptible

    command {
        gunzip -c ${eggnog_db} > /app/eggnog-mapper-2.0.1/data/eggnog.db
        gunzip -c ${eggnog_db_diamond} > /app/eggnog-mapper-2.0.1/data/eggnog_proteins.dmnd

        python /app/eggnog-mapper-2.0.1/emapper.py --cpu ${eggnog_num_cores} -i ${gene_catalogue} --output nr-eggnog --output_dir . -m diamond -d none --tax_scope auto --go_evidence non-electronic --target_orthologs all --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0 --translate --override

    }

    output {
        File eggnog_annotations = "nr-eggnog.emapper.annotations"
        File eggnog_seed_orthologs = "nr-eggnog.emapper.seed_orthologs"
    }

    runtime {
        docker: "gcr.io/osullivan-lab/eggnog-mapper:v2.0.1"
        cpu: eggnog_num_cores
        memory: eggnog_memory_gb + "GB"
        preemptible: num_preemptible
        maxRetries: num_preemptible + 1
        bootDiskSizeGb: 100
        zones: "us-central1-a us-central1-b us-central1-c us-central1-f"
        disks: "local-disk 200 HDD"
    }
}
