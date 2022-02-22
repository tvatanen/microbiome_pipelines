workflow strainphlan {

	Array[File] marker_files

 	call find_clades {

 		marker_files = marker_files

 	}

	Array[String] clades = read_lines(find_clades.strainphlan_clade_list)
	scatter (clade in clades){
		
		call strainphlan_haplotypes {
			input:
			marker_files=marker_files,
			clade=clade,
		}
	}

}

task find_clades {

	Array[File] marker_files

	command <<<

        mkdir dir_marker_files
        cat ${write_lines(marker_files)} > markerfiles.txt
        while read marker_file; do
            cp $marker_file dir_marker_files
        done <markerfiles.txt

		strainphlan -s dir_marker_files/*.pkl -o . --print_clades_only > clades.txt
		cat clades.txt | cut -d: -f4 | grep s__ | sed 's/[[:space:]]//g' > clades_edited.txt

	>>>

	output {
		File strainphlan_clade_list_raw = "clades.txt"
		File strainphlan_clade_list = "clades_edited.txt"
	}

	runtime {
		docker: "gcr.io/osullivan-lab/metaphlan:v3.0.12"
		cpu: 1
  		memory: "8GB"
  		preemptible: 2
  		disks: "local-disk 500 HDD"
	}

}


task strainphlan_haplotypes {
	String clade
	Array[File] marker_files
    Int strainphlan_memory_gb


	command <<<

        mkdir dir_marker_files
        cat ${write_lines(marker_files)} > markerfiles.txt
        while read marker_file; do
            cp $marker_file dir_marker_files
        done <markerfiles.txt

		extract_markers.py --clade ${clade} -o .
		mkdir ${clade}

		strainphlan -s dir_marker_files/*.pkl -m ${clade}.fna -o ${clade} -c ${clade} --phylophlan_mode fast --nproc 4

		tar -czf ${clade}.tar.gz ${clade}
		
    >>>
	
	output {
		File strainphlan_output = "${clade}.tar.gz"
	}

	runtime {
		docker: "gcr.io/osullivan-lab/metaphlan:v3.0.12"
		cpu: 4
        memory: strainphlan_memory_gb + "GB"
  		preemptible: 2
  		disks: "local-disk 500 HDD"
	}

}
