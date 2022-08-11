#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Defining the reference-generator workflow
workflow {
	
	ch_mhc_count = Channel.of( 0..params.mhc_allele_count )
	
	ch_kir_count = Channel.of( 0..params.kir_allele_count )
	
	ch_mhc_prot = Channel.of( 0..params.mhc_protein_count )
	
	ch_kir_prot = Channel.of( 0..params.kir_protein_count )


	PULL_IPD_MHC (
		ch_mhc_count
	)
	
	CONCAT_MHC (
		PULL_IPD_MHC.out
	)
	
	PULL_IPD_KIR (
		ch_kir_count
	)
	
	CONCAT_KIR (
		PULL_IPD_KIR.out
	)
	
	PULL_MHC_PROTEINS (
		ch_mhc_prot
	)
	
	CONCAT_MHC_PROTEINS (
		PULL_MHC_PROTEINS.out
	)
	
	PULL_KIR_PROTEINS (
		ch_kir_prot
	)
	
	CONCAT_KIR_PROTEINS (
		PULL_KIR_PROTEINS.out
	)

	CLEAN_IPD (
		CONCAT_MHC.out
		.flatten()
		.map{ file -> tuple(file.getSimpleName(), file) }
		.mix( 
			
			CONCAT_KIR.out
			.flatten()
			.map{ file -> tuple(file.getSimpleName(), file) }
		
		)
	)

	IWES_TRIMMING (
		CLEAN_IPD.out
	)

	MISEQ_TRIMMING (
		CLEAN_IPD.out
	)

	ALLELE_SORTING (
		MISEQ_TRIMMING.out
	)

	ALLELE_GROUP_NAMING (
		ALLELE_SORTING.out
	)

}


// Defining each process that will occur while generating new IPD references
process PULL_IPD_MHC {

	// This process pulls the current full roster non-human primate MHC alleles, as listed
	// in the latest Immuno Polymorphism Database release. Species-specific databases will
	// also be downloaded for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus mac-
	// aque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nem-
	// estrina, a.k.a. Mame)
	
	tag "${ipd_num}"
	
	time '1 minute'
	errorStrategy 'retry'
	maxRetries 4
	
	publishDir params.mhc_temp, pattern: '*.gbk', mode: 'move'
	
	when:
	params.pull_mhc == true
	
	input:
	val(ipd_num)

	output:
	tuple val(ipd_num), path("*.gbk")

	script:
	"""
	
	download_ipd-mhc_sequences.py ${ipd_num}
	
	"""

}


process CONCAT_MHC {
	
	when:
	ipd_num == params.mhc_allele_count
	
	input:
	tuple val(ipd_num), path(gbk)
	
	output:
	path("*.gbk")
	
	shell:
	date = new java.util.Date().format('yyyy-MM-dd')
	
	"""
	cat ${params.mhc_temp}/ipd-mhc-mafa*.gbk > ipd-mhc-mafa-${date}.gbk
	cat ${params.mhc_temp}/ipd-mhc-mamu*.gbk > ipd-mhc-mamu-${date}.gbk
	cat ${params.mhc_temp}/ipd-mhc-mane*.gbk > ipd-mhc-mane-${date}.gbk
	cat ${params.mhc_temp}/ipd-mhc-nhp*.gbk > ipd-mhc-nhp-${date}.gbk
	
	rm -rf ${params.mhc_temp}
	"""
	
}


process PULL_IPD_KIR {

	// This process pulls the current full roster non-human primate KIR alleles, as listed
	// in the latest Immuno Polymorphism Database release. Species-specific databases will
	// also be downloaded for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus mac-
	// aque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nem-
	// estrina, a.k.a. Mame)
	
	tag "${ipd_num}"
	
	time '1 minute'
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	params.pull_kir == true
	
	input:
	val(ipd_num)

	output:
	tuple val(ipd_num), path("*.gbk")

	script:
	"""

	download_ipd-kir_sequences.py ${ipd_num}

	"""

}


process CONCAT_KIR {
	
	tag "${ipd_num}"
	publishDir params.results, mode: 'move'
	
	when:
	ipd_num == params.kir_allele_count
	
	input:
	tuple val(ipd_num), path(gbk)
	
	output:
	path("*.gbk")
	
	script:
	date = new java.util.Date().format( 'yyyy-MM-dd')
	
	"""
	cat ipd-kir-mafa*.gbk > ipd-kir-mafa-${date}.gbk
	cat ipd-kir-mamu*.gbk > ipd-kir-mamu-${date}.gbk
	cat ipd-kir-mane*.gbk > ipd-kir-mane-${date}.gbk
	cat ipd-kir-nhp*.gbk > ipd-kir-nhp-${date}.gbk
	
	rm -rf ${params.kir_temp}
	"""
	
}


process PULL_MHC_PROTEINS {

	// This process pulls the current full roster non-human MHC proteins, as listed in
	// the latest Immuno Polymorphism Database release.
	
	tag "${ipd_num}"
	
	time '1 minute'
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	params.pull_mhc_proteins == true
	
	input:
	val(ipd_num)

	output:
	tuple val(ipd_num), path("*.fasta")

	script:
	"""
	download_ipd-mhc_proteins.py ${ipd_num}
	"""

}


process CONCAT_MHC_PROTEINS {
	
	tag "${ipd_num}"
	publishDir params.results, mode: 'move'
	
	when:
	ipd_num == params.mhc_protein_count
	
	input:
	tuple val(ipd_num), path(fasta)
	
	output:
	path("*.fasta")
	
	script:
	date = new java.util.Date().format( 'yyyy-MM-dd')
	
	"""
	cat ipd-mhc-mafa-prot*.gbk > ipd-mhc-mafa-prot-${date}.gbk
	cat ipd-mhc-mamu-prot*.gbk > ipd-mhc-mamu-prot-${date}.gbk
	cat ipd-mhc-mane-prot*.gbk > ipd-mhc-mane-prot-${date}.gbk
	cat ipd-mhc-nhp-prot*.gbk > ipd-mhc-nhp-prot-${date}.gbk
	
	rm -rf ${params.mhc_prot_temp}
	"""
	
}


process PULL_KIR_PROTEINS {

	// This process pulls the current full roster non-human KIR proteins, as listed
	// in the latest Immuno Polymorphism Database release.
	
	tag "${ipd_num}"
	
	time '1 minute'
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	params.pull_kir_proteins == true
	
	input:
	val(ipd_num)

	output:
	tuple val(ipd_num), path("*.fasta")

	script:
	"""

	download_ipd-kir_proteins.py ${ipd_num}

	"""

}


process CONCAT_KIR_PROTEINS {
	
	tag "${ipd_num}"
	publishDir params.results, mode: 'move'
	
	when:
	ipd_num == params.kir_protein_count
	
	input:
	tuple val(ipd_num), path(fasta)
	
	output:
	path("*.fasta")
	
	script:
	date = new java.util.Date().format( 'yyyy-MM-dd')
	
	"""
	cat ipd-kir-mafa-prot*.gbk > ipd-kir-mafa-prot-${date}.gbk
	cat ipd-kir-mamu-prot*.gbk > ipd-kir-mamu-prot-${date}.gbk
	cat ipd-kir-mane-prot*.gbk > ipd-kir-mane-prot-${date}.gbk
	cat ipd-kir-nhp-prot*.gbk > ipd-kir-nhp-prot-${date}.gbk
	
	rm -rf ${params.kir_prot_temp}
	"""
	
}


process CLEAN_IPD {

	// This process removes X's that are in amongst the bases in each sequence, and b) removes
	// sequences that are less than 100 base pairs long.

	tag "${tag}"
	publishDir params.results, mode: 'copy'

	input:
	tuple val(name), path(gbk)

	output:
	tuple val(animal_name), val(locus_name), path("*_cleaned.gbk")

	script:

	animal_name = name.substring(8,12)
	locus_name = name.substring(4,7)
	tag = name.substring(4,12)

	"""

	ipd_genbank_cleaner.py ${animal_name} ${gbk}

	"""
}


process IWES_TRIMMING {

	// This process creates databases from IPD sequences that can be used as references when
	// genotyping from immunoWES data. This means preferring genomic DNA sequences, when avail-
	// able, and falling back to exon 2 sequences when that is not an option. Trimming the data-
	// bases to exon 2 will use the same strategy I used when making miSeq amplicon trimmed data-
	// bases.

	tag "${animal_name}"
	publishDir params.results, mode: 'move'
	
	when:
	locus_name == "mhc" && animal_name == "mamu" && params.trim_for_iwes == true

	input:
	tuple val(animal_name), val(locus_name), path(gbk)

	output:
	path("*")

	script:
	"""

	trim_to_immunowes.py ${animal_name} ${gbk} ${params.iwes_exemplar}

	"""

}

process MISEQ_TRIMMING {

	// This process removes primers from IPD sequences and then deduplicates identical sequences,
	// so that groups of identical sequences can be used when genotyping.

	tag "${animal_name}"
	
	when:
	locus_name == "mhc" && animal_name == "mamu" && params.trim_for_miseq == true

	input:
	tuple val(animal_name), val(locus_name), path(gbk)

	output:
	tuple val(animal_name), path("*.miseq.trimmed.deduplicated.fasta")

	script:
	"""

	trim_to_miseq_amplicon.py ${animal_name} ${gbk} ${params.miseq_exemplar}

	"""

}


process ALLELE_SORTING {
	
	// This process sorts any lists of MHC alleles so their allele group can
	// be classified correctly
	
	tag "${animal_name}"
	
	input:
	tuple val(animal_name), path(fasta)
	
	output:
	tuple val(animal_name), path("*.sorted.fasta")
	
	script:
	"""
	
	allele_group_sorting.R ${fasta}
	
	"""
	
}


process ALLELE_GROUP_NAMING {
	
	// This process classifies allele "groups" for instances where a reference allele
	// sequence matches with numerous alleles
	
	tag "${animal_name}"
	publishDir params.results, mode: 'move'
	
	input:
	tuple val(animal_name), path(fasta)
	
	output:
	path("*")
	
	script:
	"""
	
	allele_group_naming.R ${fasta}
	
	"""
	
}
