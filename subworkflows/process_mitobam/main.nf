// mitobam subworkflow: library MT BAM + cell barcodes -> single-cell tables

nextflow.enable.dsl = 2
include { SPLIT_BATCHES } from "./modules/split_batches.nf"
include { EXTRACT_FASTA } from "./modules/extract_fasta.nf"
include { MAEGATK } from "./modules/maegatk.nf"
include { COLLAPSE_MAEGATK } from "./modules/maegatk.nf"
include { MGATK } from "./modules/mgatk.nf"
include { COLLAPSE_MGATK } from "./modules/mgatk.nf"
include { REDEEM_CONSENSUS } from "./modules/redeem.nf"
include { FINALIZE_REDEEM } from "./modules/redeem.nf"

// 

def pairBatches(batches) {
    // [sample, batch BAMs, batch barcode files] -> one [sample, batch BAM, batch barcodes] per batch
    return batches.flatMap { sample, bams, cbs ->
        def b = [bams].flatten().sort { it.name }
        def c = [cbs].flatten().sort { it.name }
        [b, c].transpose().collect { x -> [sample, x[0], x[1]] }
    }
}

//----------------------------------------------------------------------------//
// mitobam subworkflow
//----------------------------------------------------------------------------//

workflow process_mitobam {
     
    take:
        ch_mitobam
        genome

    main:

        EXTRACT_FASTA(params.string_MT, genome)

        // One pass: the library BAM becomes one BAM per batch of cells (--batch_size)
        SPLIT_BATCHES(ch_mitobam)
        ch_batches = pairBatches(SPLIT_BATCHES.out.batches)

        if (params.assay == "MAESTER") {

            // maegatk: per batch, cell by cell (UMI consensus with fgbio, re-mapping, allele counts)
            MAEGATK(ch_batches, EXTRACT_FASTA.out.fasta)
            COLLAPSE_MAEGATK(MAEGATK.out.tables.groupTuple(by:0))
            ch_output = COLLAPSE_MAEGATK.out.ch_output

        } else if (params.assay == "mtscATAC") {

            // mgatk bcall: per batch, cell by cell (mgatk's per-sample function)
            ch_ref = EXTRACT_FASTA.out.fasta.map { it -> tuple(it[0], it[5]) }
            MGATK(ch_batches, ch_ref)
            COLLAPSE_MGATK(MGATK.out.tables.groupTuple(by:0), ch_ref)
            ch_output = COLLAPSE_MGATK.out.ch_output

        } else if (params.assay == "ReDeeM") {

            // RedeemV: per batch, molecule consensus cell by cell, then strand-bias filtering of the library
            ch_ref = EXTRACT_FASTA.out.fasta.map { it -> tuple(it[0], it[5]) }
            REDEEM_CONSENSUS(ch_batches, ch_ref)
            FINALIZE_REDEEM(REDEEM_CONSENSUS.out.tables.groupTuple(by:0).map { s, g, c -> tuple(s, g.flatten(), c.flatten()) })
            ch_output = FINALIZE_REDEEM.out.ch_output

        } else {

            error "process_mitobam: no tables workflow for assay ${params.assay}"

        }

    emit:

        tables = ch_output

}

//----------------------------------------------------------------------------//
