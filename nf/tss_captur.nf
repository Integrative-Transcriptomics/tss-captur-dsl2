#!/usr/bin/env nextflow

// Info message
log.info """\

            TSS-CAPTUR
 ===================================
 MasterTable : ${params.masterTable}
 Genomes path: ${params.genomesPath}
 GFFs path   : ${params.gffPath}
 Output path : ${params.outputPath}
 """
 
include { DATAPREPARATION } from './subworkflows/dataprep'
include { MEME } from './modules/meme'
include { CLASSIFICATION } from './subworkflows/classification'
include { TERMINATORPREDICTION } from './subworkflows/terminatorpred'
include { RNAFOLD } from './modules/rnafold'
include { CREATEREPORT } from './modules/report'

workflow {
    genomesExt = Channel.fromPath(["fa", "fna", "fasta", "frn", "faa", "ffn"].collect{ "${params.inputGenomes}*.${it}" })

    DATAPREPARATION(params.masterTable, params.genomesPath, params.gffPath, params.outputPath)
    MEME(DATAPREPARATION.out.promoters, params.outputPath)
    CLASSIFICATION(DATAPREPARATION.out.queries, DATAPREPARATION.out.blastFiltered.collect(), params.outputPath)
    TERMINATORPREDICTION(genomesExt,
                        projectDir, 
                        CLASSIFICATION.out.crdFiles, 
                        DATAPREPARATION.out.summaryTranscripts, 
                        params.outputPath)
    RNAFOLD(TERMINATORPREDICTION.out.allocation, params.genomesPath, params.outputPath)
    CREATEREPORT(RNAFOLD.out.outputFigures.collect(), MEME.out.motifResult.collect(), params.outputPath)
}