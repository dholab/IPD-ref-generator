#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Defining the reference-generator workflow
workflow {

	// PULL_IPD_MHC ()
	
	// PULL_IPD_KIR ()
	
	PULL_MHC_PROTEINS ()
	
	PULL_KIR_PROTEINS ()

	// CLEAN_IPD (
	// 	PULL_IPD_MHC.out
	// 	.flatten()
	// 	.map{ file -> tuple(file.getSimpleName(), file) }
	// 	.mix( 
	// 		
	// 		PULL_IPD_KIR.out
	// 		.flatten()
	// 		.map{ file -> tuple(file.getSimpleName(), file) }
	// 	
	// 	)
	// )

// 	IWES_TRIMMING (
// 		CLEAN_IPD.out
// 	)
//
// 	MISEQ_TRIMMING (
// 		CLEAN_IPD.out
// 	)
// 
// 	ALLELE_SORTING (
// 		MISEQ_TRIMMING.out
// 	)
// 
// 	ALLELE_GROUP_NAMING (
// 		ALLELE_SORTING.out
// 	)

}


// Defining each process that will occur while generating new IPD references
// process PULL_IPD_MHC {
// 
// 	// This process pulls the current full roster non-human primate MHC alleles, as listed
// 	// in the latest Immuno Polymorphism Database release. Species-specific databases will
// 	// also be downloaded for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus mac-
// 	// aque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nem-
// 	// estrina, a.k.a. Mame)
// 
// 	time '4hours'
// 	
// 	when:
// 	params.pull_mhc == true
// 
// 	output:
// 	path("*.gbk")
// 
// 	script:
// 	"""
// 
// 	download_ipd-mhc_sequences.py ${params.mhc_allele_count}
// 
// 	"""
// 
// }


// process PULL_IPD_KIR {
// 
// 	// This process pulls the current full roster non-human primate KIR alleles, as listed
// 	// in the latest Immuno Polymorphism Database release. Species-specific databases will
// 	// also be downloaded for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus mac-
// 	// aque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nem-
// 	// estrina, a.k.a. Mame)
// 
// 	time '4hours'
// 	
// 	when:
// 	params.pull_kir == true
// 
// 	output:
// 	path("*.gbk")
// 
// 	script:
// 	"""
// 
// 	download_ipd-kir_sequences.py ${params.kir_allele_count}
// 
// 	"""
// 
// }


process PULL_MHC_PROTEINS {

	// This process pulls the current full roster non-human MHC proteins, as listed in
	// the latest Immuno Polymorphism Database release.
	
	publishDir params.results, mode: 'move'

	time '4hours'
	
	when:
	params.pull_mhc_proteins == true

	output:
	path("*.gbk")

	script:
	"""

	download_ipd-mhc_proteins.py ${params.mhc_protein_count}

	"""

}


process PULL_KIR_PROTEINS {

	// This process pulls the current full roster non-human KIR proteins, as listed
	// in the latest Immuno Polymorphism Database release.
	
	publishDir params.results, mode: 'move'

	time '4hours'
	
	when:
	params.pull_kir_proteins == true

	output:
	path("*.gbk")

	script:
	"""

	download_ipd-kir_proteins.py ${params.kir_protein_count}

	"""

}


// process CLEAN_IPD {
// 
// 	// This process removes X's that are in amongst the bases in each sequence, and b) removes
// 	// sequences that are less than 100 base pairs long.
// 
// 	tag "${animal_name}"
// 	publishDir params.results, mode: 'copy'
// 
// 	input:
// 	tuple val(name), path(gbk)
// 
// 	output:
// 	tuple val(animal_name), path("*._cleaned.gbk")
// 
// 	script:
// 
// 	animal_name = name.substring(8,12)
// 
// 	"""
// 
// 	ipd_genbank_cleaner.py ${animal_name} ${gbk}
// 
// 	"""
// }

// process IWES_TRIMMING {
// 
// 	// This process creates databases from IPD sequences that can be used as references when
// 	// genotyping from immunoWES data. This means preferring genomic DNA sequences, when avail-
// 	// able, and falling back to exon 2 sequences when that is not an option. Trimming the data-
// 	// bases to exon 2 will use the same strategy I used when making miSeq amplicon trimmed data-
// 	// bases.
// 
// 	tag "${animal_name}"
// 	publishDir params.results, mode: 'move'
// 	
// 	when:
// 	prot != "prot"
// 
// 	input:
// 	tuple val(animal_name), val(prot), path(gbk)
// 
// 	output:
// 	path("*")
// 
// 	script:
// 	"""
// 
// 	trim_to_immunowes.py ${animal_name} ${gbk} ${params.iwes_exemplar}
// 
// 	"""
// 
// }
// 
// process MISEQ_TRIMMING {
// 
// 	// This process removes primers from IPD sequences and then deduplicates identical sequences,
// 	// so that groups of identical sequences can be used when genotyping.
// 
// 	tag "${animal_name}"
// 	
// 	when:
// 	prot != "prot"
// 
// 	input:
// 	tuple val(animal_name), val(prot), path(gbk)
// 
// 	output:
// 	tuple val(animal_name), path("*.miseq.trimmed.deduplicated.fasta")
// 
// 	script:
// 	"""
// 
// 	trim_to_miseq_amplicon.py ${animal_name} ${gbk} ${params.miseq_exemplar}
// 
// 	"""
// 
// }
// 
// 
// process ALLELE_SORTING {
// 	
// 	tag ${animal_name}
// 	
// 	input:
// 	tuple val(animal_name), path(fasta)
// 	
// 	output:
// 	tuple val(animal_name), path("*.sorted.fasta")
// 	
// 	script:
// 	"""
// 	
// 	allele_group_sorting.R ${fasta}
// 	
// 	"""
// 	
// }
// 
// 
// process ALLELE_GROUP_NAMING {
// 	
// 	tag ${animal_name}
// 	publishDir params.results, mode: 'move'
// 	
// 	input:
// 	tuple val(animal_name), path(fasta)
// 	
// 	output:
// 	path("*")
// 	
// 	script:
// 	"""
// 	
// 	allele_group_naming.R ${fasta}
// 	
// 	"""
// 	
// }
