#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



workflow {
	
	PULL_IPD ()
	
	CLEAN_IPD (
		PULL_IPD.out.flatten().map{ file -> tuple(file.getSimpleName(), file) }
	)
	
	IWES_TRIMMING (
		CLEAN_IPD.out
	)
	
	MISEQ_TRIMMING (
		CLEAN_IPD.out
	)
	
}



process PULL_IPD {
	
	output:
	path("*.gbk")
	
	script:
	"""
	
	download_ipd_sequences.py ${params.allele_count}
	
	"""
	
}

process CLEAN_IPD {
	
	tag "${animal_name}"
	publishDir params.results, mode: 'copy'
	
	input:
	tuple val(name), path(gbk)
	
	output:
	tuple val(animal_name), path("*._cleaned.gbk")
	
	script:
	
	animal_name = name.substring(9,12)
	
	"""
	
	ipd_genbank_cleaner.py ${name} ${gbk}
	
	"""
}

process IWES_TRIMMING {
	
	tag "${animal_name}"
	publishDir params.results, mode: 'move'
	
	input:
	tuple val(animal_name), path(gbk)
	
	output:
	path("*")
	
	script:
	"""
	
	trim_to_immunowes.py ${name} ${gbk}
	
	"""
	
}

process MISEQ_TRIMMING {
	
	tag "${animal_name}"
	publishDir params.results, mode: 'move'
	
	input:
	tuple val(animal_name), path(gbk)
	
	output:
	path("*")
	
	script:
	"""
	
	trim_to_miseq_amplicon.py ${name} ${gbk}
	
	"""
	
}
