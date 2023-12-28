#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { badread_simulate as simulate_nuc_reads }   from './modules/simulate_reads.nf'
include { badread_simulate as simulate_mit_reads }   from './modules/simulate_reads.nf'
include { simulate_contaminant_reads }           from './modules/simulate_reads.nf'
include { downsample_simulated_reads }           from './modules/simulate_reads.nf'
include { downsample_contaminant_reads }         from './modules/simulate_reads.nf'
include { introduce_contaminants }               from './modules/simulate_reads.nf'
include { fastp }                                from './modules/simulate_reads.nf'
include { minimap2_align }                       from './modules/simulate_reads.nf'
include { qualimap_bamqc }                       from './modules/simulate_reads.nf'
include { samtools_stats }                       from './modules/simulate_reads.nf'
include { combine_alignment_qc }                 from './modules/simulate_reads.nf'
include { combine_reads }                 from './modules/simulate_reads.nf'
include { combine_nuc_mito }                 from './modules/simulate_reads.nf'


workflow {
    // Print parameters
    println "Nuclear assembly path: ${params.nuc_assembly}"
    println "Mitochondrial assembly path: ${params.mito_assembly}"
    println "Nuclear depth: ${params.nuc_depth}"
    println "Mitochondrial depth: ${params.mito_depth}"
    nuc_fasta_ch = Channel.fromPath(params.nuc_assembly).map{ it -> [it.baseName, it] }.unique{ it -> it[0] }
    mit_fasta_ch = Channel.fromPath(params.mito_assembly).map{ it -> [it.baseName, it] }.unique{ it -> it[0] }
    nuc_depth_ch = Channel.of(params.nuc_depth)
    mito_dept_ch = Channel.of(params.mito_depth)

    ch_replicates = Channel.fromList([1..params.replicates][0])

    nuclear_info_ch = nuc_fasta_ch.combine(nuc_depth_ch).combine(ch_replicates)
    mitocho_info_ch = mit_fasta_ch.combine(mito_dept_ch).combine(ch_replicates)
    
/*    nuclear_info_ch.view()
    mitocho_info_ch.view() */

    
	nuclear_sim_ch = simulate_nuc_reads(nuclear_info_ch)
	mitochon_sim_ch = simulate_mit_reads(mitocho_info_ch)
    
    nucl_reads_ch = nuclear_sim_ch.reads
    mito_reads_ch = mitochon_sim_ch.reads

    combine_reads(nucl_reads_ch, mito_reads_ch)
    combine_nuc_mito(nuc_fasta_ch,mit_fasta_ch,nucl_reads_ch.map{ it[1] })
    
    combined_reads_ch = combine_reads.out.untrimmed_reads
    combined_ref_ch = combine_nuc_mito.out.combined_ref_fasta


    combined_ref_ch.cross(combined_reads_ch).map{ it -> [it[1][0], it[1][1], it[1][2], it[0][2]] }.view()
    fastp(combined_reads_ch)
    minimap2_align(combined_ref_ch.cross(combined_reads_ch).map{ it -> [it[1][0], it[1][1], it[1][2], it[0][2]] })
    samtools_stats(minimap2_align.out)
    qualimap_bamqc(minimap2_align.out)
    combine_alignment_qc(qualimap_bamqc.out.alignment_qc.join(samtools_stats.out.stats_summary_csv, by: [0, 1]))

}
