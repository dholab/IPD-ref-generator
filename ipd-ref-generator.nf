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
		.splitCsv( header: true )
		.map { row -> tuple(row.accession, row.formal_name, row.informal_name, row.ipd_id, row.animal_id) }
	
	
	PULL_SPEC_SAMPLES (
		ch_added_seqs
	)
	
	CONCAT_SPEC_SAMPLES (
		PULL_SPEC_SAMPLES.out.collect()
	)
	
	PULL_IPD_HLA (
		ch_hla_count
	)
	
	CONCAT_HLA (
		PULL_IPD_HLA.out
	)

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

	IWES_TRIMMING (
		CLEAN_ALLELES.out
	)

	MISEQ_TRIMMING (
		CLEAN_ALLELES.out
		.mix(
			
			CONCAT_HLA.out
			.map{ file -> tuple("hum", "hla", file) }
			
		)
	)

	ALLELE_SORTING (
		MISEQ_TRIMMING.out
	)

	ALLELE_GROUP_NAMING (
		ALLELE_SORTING.out
	)

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
params.gbk_additions_temp = params.results + "/" + "gbk_add_tmp"

// iWES results
params.iwes_results = params.results + "/" + "iwes_databases"

// miseq results
params.miseq_results = params.results + "/" + "miseq_databases"


// Defining each process that will occur while generating new IPD references
process PULL_SPEC_SAMPLES {
	
	time '1minute'
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	params.pull_added_seqs == true
	
	input:
	tuple val(accession), val(formal_name), val(informal_name), val(ipd_id), val(animal_id)
	
	output:
	path("*.gbk")
	
	script:
	"""
	download_additional_embl_sequences.py "${accession}" "${formal_name}"
	"""
	
}


process CONCAT_SPEC_SAMPLES {
	
	publishDir params.results, pattern: '*.gbk', mode: params.publishMode
	
	when:
	params.pull_added_seqs == true
	
	input:
	path(gbk_files)
	
	output:
	path("*.gbk")
	
	script:
	date = new java.util.Date().format('yyyy-MM-dd')
	
	"""
	cat ipd-mhc-mafa*.gbk > ipd-mhc-mafa-${date}_added.gbk
	cat ipd-mhc-mamu*.gbk > ipd-mhc-mamu-${date}_added.gbk
	cat ipd-mhc-mane*.gbk > ipd-mhc-mane-${date}_added.gbk
	cat ipd-mhc-nhp*.gbk > ipd-mhc-nhp-${date}_added.gbk
	
	find . -name "*.gbk" -size 0 -print -delete
	"""
	
}


process PULL_IPD_MHC {

	// This process pulls the current full roster non-human primate MHC alleles, as listed
	// in the latest Immuno Polymorphism Database release. Species-specific databases will
	// also be downloaded for Rhesus macaque (Macaca mulatta, a.k.a. Mamu), Cynomolgus mac-
	// aque (Macaca fascicularis a.k.a. Mafa), and Southern pig-tailed macaque (Macaca nem-
	// estrina, a.k.a. Mame)
	
	tag "${ipd_num}"
	
	time '1minute'
	errorStrategy 'retry'
	maxRetries 4
	
	publishDir params.mhc_temp, pattern: '*.gbk', mode: params.publishMode
	
	when:
	params.pull_mhc == true
	
	input:
	val(ipd_num)

	output:
	tuple val(gbk_count), path("*.gbk")

	script:
	gbk_count = new File(params.mhc_temp)
		.listFiles()
		.findAll { it.name ==~ /.*.gbk/ }
		.size() + 1
	
	"""
	
	download_ipd-mhc_sequences.py ${ipd_num}
	
	"""

}


process CONCAT_MHC {
	
	when:
	gbk_count == params.mhc_allele_count
	
	input:
	tuple val(gbk_count), path(gbk)
	
	output:
	path("*.gbk")
	
	shell:
	date = new java.util.Date().format( 'yyyy-MM-dd')
	
	'''
	
	touch ipd-mhc-mafa-!{date}.gbk
	find !{params.mhc_temp} -maxdepth 1 -type f -name "ipd-mhc-mafa*.gbk" > mhc_mafa_list.txt && \
	for i in $(cat mhc_mafa_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_temp}/$f >> ipd-mhc-mafa-!{date}.gbk
		if test -f !{params.results}/ipd-mhc-mafa-!{date}_added.gbk; then
			cat !{params.results}/ipd-mhc-mafa-!{date}_added.gbk >> ipd-mhc-mafa-!{date}.gbk
			rm !{params.results}/ipd-mhc-mafa-!{date}_added.gbk
		fi
	done
	
	touch ipd-mhc-mamu-!{date}.gbk
	find !{params.mhc_temp} -maxdepth 1 -type f -name "ipd-mhc-mamu*.gbk" > mhc_mamu_list.txt && \
	for i in $(cat mhc_mamu_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_temp}/$f >> ipd-mhc-mamu-!{date}.gbk
		if test -f !{params.results}/ipd-mhc-mamu-!{date}_added.gbk; then
			cat !{params.results}/ipd-mhc-mamu-!{date}_added.gbk >> ipd-mhc-mamu-!{date}.gbk
			rm !{params.results}/ipd-mhc-mamu-!{date}_added.gbk
		fi
	done
	
	touch ipd-mhc-mane-!{date}.gbk
	find !{params.mhc_temp} -maxdepth 1 -type f -name "ipd-mhc-mane*.gbk" > mhc_mane_list.txt && \
	for i in $(cat mhc_mane_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_temp}/$f >> ipd-mhc-mane-!{date}.gbk
		if test -f !{params.results}/ipd-mhc-mane-!{date}_added.gbk; then
			cat !{params.results}/ipd-mhc-mane-!{date}_added.gbk >> ipd-mhc-mane-!{date}.gbk
			rm !{params.results}/ipd-mhc-mane-!{date}_added.gbk
		fi
	done
	
	touch ipd-mhc-nhp-!{date}.gbk
	find !{params.mhc_temp} -maxdepth 1 -type f -name "ipd-mhc-nhp*.gbk" > mhc_nhp_list.txt && \
	for i in $(cat mhc_nhp_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_temp}/$f >> ipd-mhc-nhp-!{date}.gbk
		if test -f !{params.results}/ipd-mhc-nhp*!{date}_added.gbk; then
			cat !{params.results}/ipd-mhc-nhp*!{date}_added.gbk >> ipd-mhc-nhp-!{date}.gbk
			rm !{params.results}/ipd-mhc-nhp*!{date}_added.gbk
		fi
	done
	
	rm -rf !{params.mhc_temp}
	find . -name "*.gbk" -size 0 -print -delete
	
	'''
	
}


process PULL_IPD_HLA {

	// This process pulls the current full roster HLA alleles, as listed
	// in the latest Immuno Polymorphism Database release.
	
	tag "${ipd_num}"
	
	cpus 1
	time '2minutes'
	errorStrategy 'retry'
	maxRetries 4
	
	publishDir params.hla_temp, pattern: '*.gbk', mode: params.publishMode
	
	when:
	params.pull_hla == true
	
	input:
	val(ipd_num)

	output:
	tuple val(gbk_count), path("*.gbk")

	script:
	gbk_count = new File(params.hla_temp)
		.listFiles()
		.findAll { it.name ==~ /.*.gbk/ }
		.size() + 1
	
	"""
	
	download_ipd-hla_sequences.py ${ipd_num}
	
	"""

}


process CONCAT_HLA {
	
	publishDir params.results, pattern: '*.gbk', mode: params.publishMode
	
	when:
	gbk_count == params.hla_allele_count
	
	input:
	tuple val(gbk_count), path(gbk)
	
	output:
	path("*.gbk")
	
	shell:
	date = new java.util.Date().format('yyyy-MM-dd')
	
	'''
	
	touch ipd-hla-!{date}.gbk
	find !{params.hla_temp} -maxdepth 1 -type f -name "*.gbk" > gbk_list.txt && \
	for i in $(cat gbk_list.txt);
	do
		f=$(basename "$i")
		cat !{params.hla_temp}/$f >> ipd-hla-!{date}.gbk
	done && \
	rm -rf !{params.hla_temp}
	
	if test -f !{params.results}/ipd-hla-!{date}_added.gbk; then
		cat !{params.results}/ipd-hla-!{date}_added.gbk >> ipd-hla-${date}.gbk
		rm !{params.results}/ipd-hla-!{date}_added.gbk
	fi
	
	find !{params.results} -name "*.gbk" -size 0 -print -delete
	
	'''
	
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
	tuple val(gbk_count), path("*.gbk")

	script:
	gbk_count = new File(params.kir_temp)
		.listFiles()
		.findAll { it.name ==~ /.*.gbk/ }
		.size() + 1
	
	"""

	download_ipd-kir_sequences.py ${ipd_num}

	"""

}


process CONCAT_KIR {
	
	publishDir params.kir_allele_results, mode: params.publishMode
		
	when:
	gbk_count == params.kir_allele_count
	
	input:
	tuple val(gbk_count), path(gbk)
	
	output:
	path("*.gbk")
	
	shell:
	date = new java.util.Date().format( 'yyyy-MM-dd')
	
	'''
	
	touch ipd-kir-mafa-!{date}.gbk
	find !{params.kir_temp} -maxdepth 1 -type f -name "ipd-kir-mafa*.gbk" > kir_mafa_list.txt && \
	for i in $(cat kir_mafa_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_temp}/$f >> ipd-kir-mafa-!{date}.gbk
	done
	
	touch ipd-kir-mamu-!{date}.gbk
	find !{params.kir_temp} -maxdepth 1 -type f -name "ipd-kir-mamu*.gbk" > kir_mamu_list.txt && \
	for i in $(cat kir_mamu_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_temp}/$f >> ipd-kir-mamu-!{date}.gbk
	done
	
	touch ipd-kir-mane-!{date}.gbk
	find !{params.kir_temp} -maxdepth 1 -type f -name "ipd-kir-mane*.gbk" > kir_mane_list.txt && \
	for i in $(cat kir_mane_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_temp}/$f >> ipd-kir-mane-!{date}.gbk
	done
	
	touch ipd-kir-nhp-!{date}.gbk
	find !{params.kir_temp} -maxdepth 1 -type f -name "ipd-kir-nhp*.gbk" > kir_nhp_list.txt && \
	for i in $(cat kir_nhp_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_temp}/$f >> ipd-kir-nhp-!{date}.gbk
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
	tuple val(fasta_count), path("*.fasta")

	script:
	fasta_count = new File(params.mhc_prot_temp)
		.listFiles()
		.findAll { it.name ==~ /.*.fasta/ }
		.size() + 1
		
	"""
	download_ipd-mhc_proteins.py ${ipd_num}
	"""

}


process CONCAT_MHC_PROTEINS {
	
	publishDir params.mhc_prot_results, mode: params.publishMode
	
	when:
	fasta_count == params.mhc_protein_count
	
	input:
	tuple val(fasta_count), path(gbk)
	
	output:
	path("*.fasta")
	
	shell:
	date = new java.util.Date().format( 'yyyy-MM-dd')
	
	'''
	
	touch ipd-mhc-mafa-prot-!{date}.fasta
	find !{params.mhc_prot_temp} -maxdepth 1 -type f -name "ipd-mhc-mafa-prot*.fasta" > mafa_prot_list.txt && \
	for i in $(cat mafa_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_prot_temp}/$f >> ipd-mhc-mafa-prot-!{date}.fasta
	done
	
	touch ipd-mhc-mamu-prot-!{date}.fasta
	find !{params.mhc_prot_temp} -maxdepth 1 -type f -name "ipd-mhc-mamu-prot*.fasta" > mamu_prot_list.txt && \
	for i in $(cat mamu_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_prot_temp}/$f >> ipd-mhc-mamu-prot-!{date}.fasta
	done
	
	touch ipd-mhc-mane-prot-!{date}.fasta
	find !{params.mhc_prot_temp} -maxdepth 1 -type f -name "ipd-mhc-mane-prot*.fasta" > mane_prot_list.txt && \
	for i in $(cat mane_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_prot_temp}/$f >> ipd-mhc-mane-prot-!{date}.fasta
	done
	
	touch ipd-mhc-nhp-prot-!{date}.fasta
	find !{params.mhc_prot_temp} -maxdepth 1 -type f -name "ipd-mhc-nhp-prot*.fasta" > nhp_prot_list.txt && \
	for i in $(cat nhp_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.mhc_prot_temp}/$f >> ipd-mhc-nhp-prot-!{date}.fasta
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
	tuple val(ipd_num), path("*.fasta")

	script:
	fasta_count = new File(params.kir_prot_temp)
		.listFiles()
		.findAll { it.name ==~ /.*.fasta/ }
		.size() + 1
	
	"""

	download_ipd-kir_proteins.py ${ipd_num}

	"""

}


process CONCAT_KIR_PROTEINS {
	
	publishDir params.kir_prot_results, mode: params.publishMode
	
	when:
	fasta_count == params.kir_protein_count
	
	input:
	tuple val(fasta_count), path(fasta)
	
	output:
	path("*.fasta")
	
	shell:
	date = new java.util.Date().format( 'yyyy-MM-dd')
	
	'''
	
	touch ipd-kir-mafa-prot-!{date}.fasta
	find !{params.kir_prot_temp} -maxdepth 1 -type f -name "ipd-kir-mafa-prot*.fasta" > mafa_prot_list.txt && \
	for i in $(cat mafa_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_prot_temp}/$f >> ipd-kir-mafa-prot-!{date}.fasta
	done
	
	touch ipd-kir-mamu-prot-!{date}.fasta
	find !{params.kir_prot_temp} -maxdepth 1 -type f -name "ipd-kir-mamu-prot*.fasta" > mamu_prot_list.txt && \
	for i in $(cat mamu_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_prot_temp}/$f >> ipd-kir-mamu-prot-!{date}.fasta
	done
	
	touch ipd-kir-mane-prot-!{date}.fasta
	find !{params.kir_prot_temp} -maxdepth 1 -type f -name "ipd-kir-mane-prot*.fasta" > mane_prot_list.txt && \
	for i in $(cat mane_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_prot_temp}/$f >> ipd-kir-mane-prot-!{date}.fasta
	done
	
	touch ipd-kir-nhp-prot-!{date}.fasta
	find !{params.kir_prot_temp} -maxdepth 1 -type f -name "ipd-kir-nhp-prot*.fasta" > nhp_prot_list.txt && \
	for i in $(cat nhp_prot_list.txt);
	do
		f=$(basename "$i")
		cat !{params.kir_prot_temp}/$f >> ipd-kir-nhp-prot-!{date}.fasta
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

	input:
	tuple val(name), path(gbk)

	output:
	tuple val(animal_name), val(locus_name), path("*_cleaned.gbk")

	script:

	animal_name = name.substring(8,12)
	locus_name = name.substring(4,7)
	tag = name.substring(4,12)

	"""

	ipd_genbank_cleaner.py ${animal_name} ${locus_name} ${gbk}

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
	params.trim_for_iwes == true && locus_name == "mhc" && (animal_name == "mamu" || animal_name == "mafa")

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
	publishDir params.miseq_results, pattern: "*hla*.fasta", mode: 'copy'
	
	when:
	params.trim_for_miseq == true && locus_name != "kir" && animal_name != "nhp"

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
