#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Defining the reference-generator workflow
workflow {
	
	ch_hla_count = Channel
		.of( 1..params.hla_allele_count )
		.toSortedList()
		.flatten()
	
	ch_mhc_count = Channel
		.of( 1..params.mhc_allele_count )
		.toSortedList()
		.flatten()
	
	// ch_mhc_count = Channel
	// 	.of( 1..params.mhc_allele_count )
	// 	.buffer( size: 100, remainder: true)
	
	ch_kir_count = Channel
		.of( 1..params.kir_allele_count )
		.toSortedList()
		.flatten()
	
	ch_mhc_prot = Channel
		.of( 1..params.mhc_protein_count )
		.toSortedList()
		.flatten()
	
	ch_kir_prot = Channel
		.of( 1..params.kir_protein_count )
		.toSortedList()
		.flatten()
		
	ch_added_seqs = Channel
		.fromPath( params.added_seqs )
	
	
	PULL_SPEC_SAMPLES (
		ch_added_seqs,
		ch_added_seqs
			.splitCsv()
			.count()
	)
	
	// CONCAT_SPEC_SAMPLES (
	// 	PULL_SPEC_SAMPLES.out.collect()
	// )

	PULL_IPD_MHC ()
	
	// CONCAT_MHC (
	// 	PULL_IPD_MHC.out
	// 		.flatten()
	// 		.map { it -> it.getName() }
	// 		.collectFile(name: 'allele_file_list.txt', newLine: true)
	// )

	CONCAT_MHC (
		PULL_IPD_MHC.out
			.mix (
				PULL_SPEC_SAMPLES.out
			)
			.collect()
		
	)
	
	PULL_IPD_HLA (
		ch_hla_count
	)
	
	CONCAT_HLA (
		PULL_IPD_HLA.out
	)
	
	COMPRESS_HLA (
		CONCAT_HLA.out
	)
	
	PULL_IPD_KIR (
		ch_kir_count
	)
	
	CONCAT_KIR (
		PULL_IPD_KIR.out
			.flatten()
			.map { it -> it.getName() }
			.collectFile(name: 'allele_file_list.txt', newLine: true)
	)
	
	PULL_MHC_PROTEINS (
		ch_mhc_prot
	)
	
	CONCAT_MHC_PROTEINS (
		PULL_MHC_PROTEINS.out
			.flatten()
			.map { it -> it.getName() }
			.collectFile(name: 'allele_file_list.txt', newLine: true)
	)
	
	PULL_KIR_PROTEINS (
		ch_kir_prot
	)
	
	CONCAT_KIR_PROTEINS (
		PULL_KIR_PROTEINS.out
			.flatten()
			.map { it -> it.getName() }
			.collectFile(name: 'allele_file_list.txt', newLine: true)
	)

	// CLEAN_ALLELES (
	// 	CONCAT_MHC.out
	// 	.flatten()
	// 	.map{ file -> tuple(file.getSimpleName(), file) }
	// 	.mix( 
			
	// 		CONCAT_KIR.out
	// 		.flatten()
	// 		.map{ file -> tuple(file.getSimpleName(), file) }
		
	// 	)
	// )

	// IWES_TRIMMING (
	// 	CLEAN_ALLELES.out
	// )

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

	// MISEQ_TRIMMING (
	// 	CLEAN_ALLELES.out
	// 	.mix(
			
	// 		CONCAT_HLA.out
	// 		.map{ file -> tuple("hum", "hla", file) }
			
	// 	)
	// )

	ALLELE_SORTING (
		MISEQ_TRIMMING.out
	)

	// ALLELE_GROUP_NAMING (
	// 	ALLELE_SORTING.out
	// )

}


// specifying whether to run in resumable mode
if( params.resumable == true ) {
	params.publishMode = 'copy'
}
else {
	params.publishMode = 'move'
}

// Defining where to place results
// hla alleles
params.hla_temp = params.results + "/" + "hla_tmp"
params.hla_results = params.results + "/" + "hla_alleles"

// mhc alleles
params.mhc_temp = params.results + "/" + "mhc_tmp"
params.mhc_allele_results = params.results + "/" + "mhc_alleles"

// kir alleles
params.kir_temp = params.results + "/" + "kir_tmp"
params.kir_allele_results = params.results + "/" + "kir_alleles"

// mhc proteins
params.mhc_prot_temp = params.results + "/" + "mhc_prot_tmp"
params.mhc_prot_results = params.results + "/" + "mhc_proteins"

// kir proteins
params.kir_prot_temp = params.results + "/" + "kir_prot_tmp"
params.kir_prot_results = params.results + "/" + "kir_proteins"

// additions from samplesheet
// params.gbk_additions_temp = params.results + "/" + "gbk_add_tmp"
params.spec_results = params.results + "/" + "samplesheet_sequences"

// Exon 2 results
params.exon2_results = params.results + "/" + "exon2_databases"

// iWES results
params.iwes_results = params.results + "/" + "iwes_databases"

// miseq results
params.miseq_results = params.results + "/" + "miseq_databases"


// Defining each process that will occur while generating new IPD references
process PULL_SPEC_SAMPLES {
	
	// publishDir params.spec_results, mode: 'copy'

	tag "${sample_count} seqs"
	
	// errorStrategy 'retry'
	// maxRetries 4
	
	when:
	params.pull_added_seqs == true
	
	input:
	path samplesheet
	val sample_count
	
	output:
	path "*.gbk"
	
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

	tag "${params.mhc_allele_count} alleles"
	
	// publishDir params.mhc_allele_results, mode: params.publishMode, overwrite: true
	
	when:
	params.pull_mhc == true
	
	output:
	path "*.gbk"

	script:
	"""
	goDownloadIPD MHC ${params.mhc_allele_count} ${params.last_release_date}
	"""

}


process CONCAT_MHC {

	publishDir params.mhc_allele_results, mode: 'copy', overwrite: true
	publishDir params.resources, pattern: '*nhp*.gbk', saveAs: "previous_mhc_database.gbk", mode: 'copy', overwrite: true
	
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
	
	tag "${ipd_num}"
	
	cpus 1
	time { 2.minute * task.attempt }
	errorStrategy 'retry'
	maxRetries 4
	
	publishDir params.hla_temp, pattern: '*.gbk', mode: params.publishMode
	
	when:
	params.pull_hla == true
	
	input:
	val(ipd_num)

	output:
	path("*.gbk")

	script:
	"""
	goDownloadIPD HLA ${ipd_num} ${params.last_release_date}
	"""

}


process CONCAT_HLA {
	
	input:
	path(gbk)
	
	output:
	path("*.gbk")
	
	shell:
	'''
	
	touch ipd-hla-!{params.date}.gbk
	find !{params.hla_temp} -maxdepth 1 -type f -name "*.gbk" > gbk_list.txt && \
	for i in $(cat gbk_list.txt);
	do
		f=$(basename "$i")
		cat !{params.hla_temp}/$f >> ipd-hla-!{params.date}.gbk
	done && \
	rm -rf !{params.hla_temp}
	
	if test -f !{params.results}/ipd-hla-!{params.date}_added.gbk; then
		cat !{params.results}/ipd-hla-!{params.date}_added.gbk >> ipd-hla-${params.date}.gbk
		rm !{params.results}/ipd-hla-!{params.date}_added.gbk
	fi
	
	find !{params.results} -name "*.gbk" -size 0 -print -delete
	
	'''
	
}


process COMPRESS_HLA {
	
	publishDir params.hla_results, pattern: '*.gbk.xz', mode: 'move'
	
	input:
	path(gbk)
	
	output:
	path(gbk)
	
	script:
	"""
	xz -9 ${gbk}
	"""
	
}


process PULL_IPD_KIR {

	// This process pulls the current full roster non-human primate KIR alleles, as listed
	// in the latest Immuno Polymorphism Database release. Species-specific databases will
	// also be downloaded for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus mac-
	// aque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nem-
	// estrina, a.k.a. Mame)
	
	tag "${ipd_num}"
	
	time '1minute'
	errorStrategy 'retry'
	maxRetries 4
	
	publishDir params.kir_temp, pattern: '*.gbk', mode: params.publishMode
	
	when:
	params.pull_kir == true
	
	input:
	val(ipd_num)

	output:
	path("*.gbk")

	script:
	"""
	goDownloadIPD KIR ${params.kir_allele_count} ${params.last_release_date}
	"""

}


process CONCAT_KIR {
	
	// publishDir params.kir_allele_results, mode: params.publishMode
	
	input:
	path(gbk)
	
	output:
	path("*.gbk")
	
	shell:
	'''
	
	touch ipd-kir-mafa-!{params.date}.gbk
	find !{params.kir_temp} -maxdepth 1 -type f -name "ipd-kir-mafa*.gbk" > kir_mafa_list.txt && \
	for i in $(cat kir_mafa_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_temp}/$f >> ipd-kir-mafa-!{params.date}.gbk
	done
	
	touch ipd-kir-mamu-!{params.date}.gbk
	find !{params.kir_temp} -maxdepth 1 -type f -name "ipd-kir-mamu*.gbk" > kir_mamu_list.txt && \
	for i in $(cat kir_mamu_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_temp}/$f >> ipd-kir-mamu-!{params.date}.gbk
	done
	
	touch ipd-kir-mane-!{params.date}.gbk
	find !{params.kir_temp} -maxdepth 1 -type f -name "ipd-kir-mane*.gbk" > kir_mane_list.txt && \
	for i in $(cat kir_mane_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_temp}/$f >> ipd-kir-mane-!{params.date}.gbk
	done
	
	touch ipd-kir-nhp-!{params.date}.gbk
	find !{params.kir_temp} -maxdepth 1 -type f -name "ipd-kir-nhp*.gbk" > kir_nhp_list.txt && \
	for i in $(cat kir_nhp_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_temp}/$f >> ipd-kir-nhp-!{params.date}.gbk
	done
	
	rm -rf !{params.kir_temp}
	find . -name "*.gbk" -size 0 -print -delete
	
	'''

	
}


process PULL_MHC_PROTEINS {

	// This process pulls the current full roster non-human MHC proteins, as listed in
	// the latest Immuno Polymorphism Database release.
	
	tag "${ipd_num}"
	
	time '1minute'
	errorStrategy 'retry'
	maxRetries 4
	
	publishDir params.mhc_prot_temp, pattern: '*.fasta', mode: params.publishMode
	
	when:
	params.pull_mhc_proteins == true
	
	input:
	val(ipd_num)

	output:
	path("*.fasta")

	script:
	"""
	
	download_ipd-mhc_proteins.py ${ipd_num}
	
	"""

}


process CONCAT_MHC_PROTEINS {
	
	publishDir params.mhc_prot_results, mode: params.publishMode
	
	input:
	path(fasta)
	
	output:
	path("*.fasta")
	
	shell:
	'''
	
	touch ipd-mhc-mafa-prot-!{params.date}.fasta
	find !{params.mhc_prot_temp} -maxdepth 1 -type f -name "ipd-mhc-mafa-prot*.fasta" > mafa_prot_list.txt && \
	for i in $(cat mafa_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_prot_temp}/$f >> ipd-mhc-mafa-prot-!{params.date}.fasta
	done
	
	touch ipd-mhc-mamu-prot-!{params.date}.fasta
	find !{params.mhc_prot_temp} -maxdepth 1 -type f -name "ipd-mhc-mamu-prot*.fasta" > mamu_prot_list.txt && \
	for i in $(cat mamu_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_prot_temp}/$f >> ipd-mhc-mamu-prot-!{params.date}.fasta
	done
	
	touch ipd-mhc-mane-prot-!{params.date}.fasta
	find !{params.mhc_prot_temp} -maxdepth 1 -type f -name "ipd-mhc-mane-prot*.fasta" > mane_prot_list.txt && \
	for i in $(cat mane_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_prot_temp}/$f >> ipd-mhc-mane-prot-!{params.date}.fasta
	done
	
	touch ipd-mhc-nhp-prot-!{params.date}.fasta
	find !{params.mhc_prot_temp} -maxdepth 1 -type f -name "ipd-mhc-nhp-prot*.fasta" > nhp_prot_list.txt && \
	for i in $(cat nhp_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_prot_temp}/$f >> ipd-mhc-nhp-prot-!{params.date}.fasta
	done
	
	rm -rf !{params.mhc_prot_temp}
	find . -name "-prot*.fasta" -size 0 -print -delete
	
	'''
	
}


process PULL_KIR_PROTEINS {

	// This process pulls the current full roster non-human KIR proteins, as listed
	// in the latest Immuno Polymorphism Database release.
	
	tag "${ipd_num}"
	
	time '1minute'
	errorStrategy 'retry'
	maxRetries 4
	
	publishDir params.kir_prot_temp, pattern: '*.fasta', mode: params.publishMode
	
	when:
	params.pull_kir_proteins == true
	
	input:
	val(ipd_num)

	output:
	path("*.fasta")

	script:
	"""

	download_ipd-kir_proteins.py ${ipd_num}

	"""

}


process CONCAT_KIR_PROTEINS {
	
	publishDir params.kir_prot_results, mode: params.publishMode
	
	input:
	path(fasta)
	
	output:
	path("*.fasta")
	
	shell:
	'''
	
	touch ipd-kir-mafa-prot-!{params.date}.fasta
	find !{params.kir_prot_temp} -maxdepth 1 -type f -name "ipd-kir-mafa-prot*.fasta" > mafa_prot_list.txt && \
	for i in $(cat mafa_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_prot_temp}/$f >> ipd-kir-mafa-prot-!{params.date}.fasta
	done
	
	touch ipd-kir-mamu-prot-!{params.date}.fasta
	find !{params.kir_prot_temp} -maxdepth 1 -type f -name "ipd-kir-mamu-prot*.fasta" > mamu_prot_list.txt && \
	for i in $(cat mamu_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_prot_temp}/$f >> ipd-kir-mamu-prot-!{params.date}.fasta
	done
	
	touch ipd-kir-mane-prot-!{params.date}.fasta
	find !{params.kir_prot_temp} -maxdepth 1 -type f -name "ipd-kir-mane-prot*.fasta" > mane_prot_list.txt && \
	for i in $(cat mane_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_prot_temp}/$f >> ipd-kir-mane-prot-!{params.date}.fasta
	done
	
	touch ipd-kir-nhp-prot-!{params.date}.fasta
	find !{params.kir_prot_temp} -maxdepth 1 -type f -name "ipd-kir-nhp-prot*.fasta" > nhp_prot_list.txt && \
	for i in $(cat nhp_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_prot_temp}/$f >> ipd-kir-nhp-prot-!{params.date}.fasta
	done
	
	rm -rf !{params.kir_prot_temp}
	find . -name "-prot*.fasta" -size 0 -print -delete
	
	'''
	
}


process CLEAN_ALLELES {

	// This process removes X's that are in amongst the bases in each sequence, and b) removes
	// sequences that are less than 100 base pairs long.

	tag "${tag}"
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
	tag = name.substring(4,12)

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
	publishDir params.exon2_results, pattern: '*exon2_deduplicated*' mode: 'move'

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

	tag "${animal_name}"
	publishDir params.iwes_results, mode: 'move'
	
	when:
	params.trim_for_iwes == true 

	input:
	tuple val(name), path(gbk)

	output:
	path "*.fasta"

	script:

	animal_name = name.substring(8,12)
	locus_name = name.substring(4,7)
	tag = name.substring(4,12)
	
	"""

	trim_to_immunowes.py ${animal_name} ${gbk} ${params.iwes_exemplar}

	"""

}


process MISEQ_TRIMMING {

	// This process removes primers from IPD sequences and then deduplicates identical sequences,
	// so that groups of identical sequences can be used when genotyping.

	tag "${animal_name}"
	publishDir params.miseq_results, pattern: "*hla*.fasta", mode: 'copy'
	
	when:
	params.trim_for_miseq == true

	input:
	tuple val(name), path(gbk)

	output:
	tuple val(animal_name), path("*.miseq.trimmed.deduplicated.fasta")

	script:

	animal_name = name.substring(8,12)
	locus_name = name.substring(4,7)
	tag = name.substring(4,12)

	"""

	trim_to_miseq_amplicon.py ${animal_name} ${gbk} ${params.miseq_exemplar}

	"""

}


process ALLELE_SORTING {
	
	// This process sorts any lists of MHC alleles so their allele group can
	// be classified correctly
	
	tag "${animal_name}"
	
	when:
	params.trim_for_miseq == true && animal_name != "nhp" && animal_name != "hum"
	
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
	publishDir params.miseq_results, mode: 'move'
	
	input:
	tuple val(animal_name), path(fasta)
	
	output:
	path("*")
	
	script:
	"""
	
	allele_group_naming.R ${fasta}
	
	"""
	
}
