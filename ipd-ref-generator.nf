#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Defining the reference-generator workflow
workflow {
	
	ch_hla_count = Channel
		.of( params.hla_allele_count )
	
	ch_mhc_count = Channel
		.of( params.mhc_allele_count )
	
	ch_kir_count = Channel
		.of( params.kir_allele_count )
	
	ch_mhc_prot = Channel
		.of( params.mhc_protein_count )
	
	ch_kir_prot = Channel
		.of( params.kir_protein_count )
		
	ch_added_seqs = Channel
		.fromPath( params.added_seqs )
	
	
	PULL_SPEC_SAMPLES (
		ch_added_seqs,
		ch_added_seqs
			.splitCsv()
			.count()
	)

	PULL_IPD_MHC (
		ch_mhc_count
	)

	CONCAT_MHC (
		PULL_IPD_MHC.out.embl.collect()
			.mix (
				PULL_SPEC_SAMPLES.out.collect()
			)
		
	)
	
	PULL_IPD_HLA (
		ch_hla_count
	)
	
	CONCAT_HLA (
		PULL_IPD_HLA.out.collect()
	)
	
	COMPRESS_HLA (
		CONCAT_HLA.out
	)
	
	PULL_IPD_KIR (
		ch_kir_count
	)
	
	CONCAT_KIR (
		PULL_IPD_KIR.out.collect()
	)
	
	PULL_MHC_PROTEINS (
		ch_mhc_prot
	)
	
	CONCAT_MHC_PROTEINS (
		PULL_MHC_PROTEINS.out.collect()
	)
	
	PULL_KIR_PROTEINS (
		ch_kir_prot
	)
	
	CONCAT_KIR_PROTEINS (
		PULL_KIR_PROTEINS.out.collect()
	)

	CLEAN_ALLELES (
		CONCAT_MHC.out
		.flatten()
		.map{ file -> tuple(file.getSimpleName(), file) }
		.mix( 
			
			CONCAT_KIR.out
			.flatten()
			.map{ file -> tuple(file.getSimpleName(), file) }
		
		)
	)

	EXON2_TRIMMING (
		CONCAT_MHC.out
			.flatten()
			.map{ file -> tuple(file.getSimpleName(), file) }
	)

	IWES_TRIMMING (
		CONCAT_MHC.out
			.flatten()
			.map{ file -> tuple(file.getSimpleName(), file) }
	)

	MISEQ_TRIMMING (
		CONCAT_MHC.out
			.flatten()
			.map{ file -> tuple(file.getSimpleName(), file) }
	)

	ALLELE_SORTING (
		MISEQ_TRIMMING.out
	)

	// ALLELE_GROUP_NAMING (
	// 	ALLELE_SORTING.out
	// )

}


// Defining where to place results
// hla alleles
params.hla_results = params.results + "/" + "hla_alleles"

// mhc alleles
params.mhc_allele_results = params.results + "/" + "mhc_alleles"

// kir alleles
params.kir_allele_results = params.results + "/" + "kir_alleles"

// mhc proteins
params.mhc_prot_results = params.results + "/" + "mhc_proteins"

// kir proteins
params.kir_prot_results = params.results + "/" + "kir_proteins"

// additions from samplesheet
params.spec_results = params.results + "/" + "samplesheet_sequences"

// Exon 2 results
params.exon2_results = params.results + "/" + "exon2_databases"

// iWES results
params.iwes_results = params.results + "/" + "iwes_databases"

// miseq results
params.miseq_results = params.results + "/" + "miseq_databases"


// Defining each process that will occur while generating new IPD references
process PULL_SPEC_SAMPLES {

	tag "${sample_count} seqs"

	cpus 1

	errorStrategy 'ignore'
	
	input:
	path samplesheet
	val sample_count
	
	output:
	path "*.gbk"
	
	when:
	params.pull_added_seqs == true
	
	script:
	"""
	download_additional_embl_sequences.py ${samplesheet} ${params.date}
	"""
	
}


process PULL_IPD_MHC {

	// This process pulls the current full roster non-human primate MHC alleles, as listed
	// in the latest Immuno Polymorphism Database release. Species-specific databases will
	// also be downloaded for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus mac-
	// aque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nem-
	// estrina, a.k.a. Mame)

	tag "${allele_count} alleles"
	// publishDir params.resources, mode: 'copy', pattern: '*.json'

	cpus params.max_shared_cpus

	errorStrategy 'ignore'

	input:
	val allele_count
	
	output:
	path "*.embl", emit: embl
	// path "*.json", emit: lookup
	
	when:
	params.pull_mhc == true

	script:
	"""
	goDownloadIPD MHC ${allele_count} ${params.last_release_date} ${params.resources}
	"""

}


process CONCAT_MHC {

	publishDir params.mhc_allele_results, mode: 'copy', overwrite: true
	// publishDir params.resources, pattern: '*nhp*.gbk', saveAs: {item -> "previous_mhc_database.gbk"}, mode: 'copy', overwrite: true
	
	cpus 1

	input:
	path gbk_files
	
	output:
	path "*.gbk"

	when:
	params.pull_mhc == true || params.pull_added_seqs == true

	script:
	"""
	collate_alleles_by_animal.py --gene MHC --previous_database "${params.resources}/previous_mhc_database.gbk"
	"""

}


process PULL_IPD_HLA {

	// This process pulls the updated roster of HLA alleles, as listed
	// in the latest Immuno Polymorphism Database release.
	
	tag "${allele_count} alleles"
	
	cpus params.max_shared_cpus

	errorStrategy 'ignore'
	
	input:
	val allele_count

	output:
	path "*.embl"
	
	when:
	params.pull_hla == true

	script:
	"""
	goDownloadIPD HLA ${allele_count} ${params.last_release_date}
	"""

}


process CONCAT_HLA {
	
	input:
	path gbk
	
	output:
	path "*.gbk"
	
	script:
	"""
	collate_alleles_by_animal.py --gene HLA --previous_database "${params.resources}/previous_hla_database.gbk"
	"""
	
}


process COMPRESS_HLA {
	
	publishDir params.hla_results, pattern: '*.zst', mode: 'copy'
	
	input:
	path gbk
	
	output:
	path "*.zst"
	
	script:
	"""
	zstd `realpath ${gbk}` -o ipd-hla-${params.date}.gbk.zst
	"""
	
}


process PULL_IPD_KIR {

	// This process pulls the current full roster non-human primate KIR alleles, as listed
	// in the latest Immuno Polymorphism Database release. Species-specific databases will
	// also be downloaded for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus mac-
	// aque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nem-
	// estrina, a.k.a. Mame)
	
	tag "${allele_count} alleles"
	
	cpus params.max_shared_cpus

	errorStrategy 'ignore'
	
	input:
	val allele_count

	output:
	path "*.embl"
	
	when:
	params.pull_kir == true

	script:
	"""
	goDownloadIPD KIR ${allele_count} ${params.last_release_date} ${params.resources}
	"""

}


process CONCAT_KIR {
	
	publishDir params.kir_allele_results, mode: 'copy'
	publishDir params.resources, pattern: '*nhp*.gbk', saveAs: "previous_kir_database.gbk", mode: 'copy', overwrite: true

	errorStrategy 'ignore'

	cpus 1
	
	input:
	path gbk
	
	output:
	path "*.gbk"
	
	script:
	"""
	collate_alleles_by_animal.py --gene KIR --previous_database "${params.resources}/previous_kir_database.gbk"
	"""

}


process PULL_MHC_PROTEINS {

	// This process pulls the current full roster non-human MHC proteins, as listed in
	// the latest Immuno Polymorphism Database release.

	tag "${protein_count} proteins"

	errorStrategy 'ignore'

	cpus params.max_shared_cpus
	
	input:
	val protein_count

	output:
	path "*.fasta"
	
	when:
	params.pull_mhc_proteins == true

	script:
	"""
	goDownloadIPD MHCPRO ${protein_count} ${params.last_release_date} ${params.resources}
	"""

}


process CONCAT_MHC_PROTEINS {
	
	publishDir params.mhc_prot_results, mode: 'copy'

	errorStrategy 'ignore'
	
	input:
	path fastas
	
	output:
	path "*.fasta"
	
	script:
	"""
	collate_proteins_by_animal.py --gene MHC
	"""
	
}


process PULL_KIR_PROTEINS {

	// This process pulls the current full roster non-human KIR proteins, as listed
	// in the latest Immuno Polymorphism Database release.
	
	tag "${protein_count} proteins"

	errorStrategy 'ignore'
	
	input:
	val protein_count

	output:
	path "*.fasta"
	
	when:
	params.pull_kir_proteins == true

	script:
	"""
	goDownloadIPD KIRPRO ${protein_count} ${params.last_release_date} ${params.resources}
	"""

}


process CONCAT_KIR_PROTEINS {
	
	publishDir params.kir_prot_results, mode: 'copy'

	errorStrategy 'ignore'
	
	input:
	path fastas
	
	output:
	path "*.fasta"
	
	script:
	"""
	collate_proteins_by_animal.py --gene KIR
	"""
	
}


process CLEAN_ALLELES {

	// This process removes X's that are in amongst the bases in each sequence, and b) removes
	// sequences that are less than 100 base pairs long.

	tag "${animal_name} ${locus_name}"
	publishDir params.mhc_allele_results, pattern: "*mhc*.gbk", mode: 'copy'
	publishDir params.kir_allele_results, pattern: "*kir*.gbk", mode: 'copy'
	publishDir params.hla_results, pattern: "*hla*.gbk", mode: 'copy'

	input:
	tuple val(name), path(gbk)

	output:
	tuple val(animal_name), val(locus_name), path("*_cleaned.gbk")

	script:
	animal_name = name.substring(8,12)
	locus_name = name.substring(4,7)
	"""
	ipd_genbank_cleaner.py --animal ${animal_name} --gene ${locus_name} --file ${gbk}
	"""
}


process EXON2_TRIMMING {

	/*
	In this process, we trim each MHC database into exon-2-only sequences, which can 
	be used for lower-resolution genotyping. These databases are useful in scenarios
	when a) high-resolution genotyping is not needed or is computationally prohibitive,
	or b) available reference databases are too incomplete to reliably perform high-
	resolution genotyping. This step creates an exon 2 database for all alleles, and
	also creates separate databases for Class I and Class II alleles, which are 
	sometimes treated separately.
	*/

	tag "${animal_name}"
	publishDir params.exon2_results, pattern: '*exon2_deduplicated*', mode: 'copy'

	input:
	tuple val(name), path(gbk)

	output:
	path "*.fasta"
	
	when:
	params.trim_for_exon2 == true 

	script:
	animal_name = name.substring(8,12)
	"""
	trim_to_exon2.py ${gbk}
	"""

}


process IWES_TRIMMING {

	// This process creates databases from IPD sequences that can be used as references when
	// genotyping from immunoWES data. This means preferring genomic DNA sequences, when avail-
	// able, and falling back to exon 2 sequences when that is not an option. Trimming the data-
	// bases to exon 2 will use the same strategy I used when making miSeq amplicon trimmed data-
	// bases.

	tag "${animal_name} ${locus_name}"
	publishDir params.iwes_results, mode: 'copy'

	input:
	tuple val(name), path(gbk)

	output:
	path "*.fasta"
	
	when:
	params.trim_for_iwes == true 

	script:
	animal_name = name.substring(8,12)
	locus_name = name.substring(4,7)
	"""
	trim_to_immunowes.py ${animal_name} ${gbk} ${params.iwes_exemplar}
	"""

}


process MISEQ_TRIMMING {

	// This process removes primers from IPD sequences and then deduplicates identical sequences,
	// so that groups of identical sequences can be used when genotyping.

	tag "${animal_name} ${locus_name}"
	publishDir params.miseq_results, pattern: "*hla*.fasta", mode: 'copy'

	input:
	tuple val(name), path(gbk)

	output:
	tuple val(animal_name), path("*.miseq.trimmed.deduplicated.fasta")
	
	when:
	params.trim_for_miseq == true

	script:
	animal_name = name.substring(8,12)
	locus_name = name.substring(4,7)
	"""
	trim_to_miseq_amplicon.py ${animal_name} ${gbk} ${params.miseq_exemplar}
	"""

}


process ALLELE_SORTING {
	
	// This process sorts any lists of MHC alleles so their allele group can
	// be classified correctly
	
	tag "${animal_name}"
	publishDir params.miseq_results, mode: 'copy'

	errorStrategy 'ignore'
	
	input:
	tuple val(animal_name), path(fasta)
	
	output:
	tuple val(animal_name), path("*.sorted.fasta")
	
	when:
	params.trim_for_miseq == true && animal_name != "nhp" && animal_name != "hum"
	
	script:
	"""
	allele_group_sorting.R ${fasta}
	"""
	
}


// process ALLELE_GROUP_NAMING {
	
// 	// This process classifies allele "groups" for instances where a reference allele
// 	// sequence matches with numerous alleles
	
// 	tag "${animal_name}"
// 	publishDir params.miseq_results, mode: 'copy'
	
// 	input:
// 	tuple val(animal_name), path(fasta)
	
// 	output:
// 	path "*"
	
// 	script:
// 	"""
// 	allele_group_naming.R ${fasta}
// 	"""
	
// }
