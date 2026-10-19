// Preprocessing subworkflows: reads (FASTQ or BAM) -> single-cell allelic tables
nextflow.enable.dsl = 2

// Subworkflows
include { prep_resources } from "../prep_resources/main"
include { prep_genome } from "../prep_resources/main"
include { get_maester_bam } from "../get_mitobam/main"
include { get_atac_bam } from "../get_mitobam/main"
include { process_mitobam } from "../process_mitobam/main"


//


def checkParams() {
    // Assay -> pre-processing tool; assays implemented in this version
    def PP_METHODS = ['MAESTER': 'maegatk', 'mtscATAC': 'mgatk', 'ReDeeM': 'redeem-v']
    def IMPLEMENTED = ['MAESTER', 'mtscATAC']
    if (!params.input)
        error "Error: --input (samplesheet CSV) is required"
    if (!PP_METHODS.containsKey(params.assay))
        error "Error: --assay must be one of ${PP_METHODS.keySet().join(', ')}. Provided: ${params.assay}"
    if (!IMPLEMENTED.contains(params.assay))
        error "Error: assay ${params.assay} is not implemented yet (available: ${IMPLEMENTED.join(', ')})"
}

// Read and validate the samplesheet. Columns: sample, fastq_dir, bam, cell_barcodes.
// Exactly one of fastq_dir / bam per row; cell_barcodes (one barcode per line) is required.
// Returns a list of [sample, fastq_dir | null, bam | null, cell_barcodes]
def readSamplesheet() {

    checkParams()
    def seen = [] as Set
    def rows = file(params.input, checkIfExists: true).splitCsv(header: true)
    if (rows.isEmpty())
        error "Samplesheet ${params.input} is empty"

    return rows.collect { row ->
        def sample = row.sample?.trim()
        def fastq = row.fastq_dir?.trim() ?: null
        def bam = row.bam?.trim() ?: null
        def barcodes = row.cell_barcodes?.trim() ?: null
        if (!sample)
            error "Samplesheet: empty 'sample' in row ${row}"
        if (!seen.add(sample))
            error "Samplesheet: duplicated sample '${sample}'"
        if ((fastq == null) == (bam == null))
            error "Samplesheet, sample '${sample}': give exactly one of 'fastq_dir' or 'bam'"
        if (barcodes == null)
            error "Samplesheet, sample '${sample}': 'cell_barcodes' is required"
        if (fastq != null && !file(fastq).isDirectory())
            error "Samplesheet, sample '${sample}': fastq_dir not found: ${fastq}"
        if (bam != null && !file(bam).exists())
            error "Samplesheet, sample '${sample}': bam not found: ${bam}"
        if (!file(barcodes).exists())
            error "Samplesheet, sample '${sample}': cell_barcodes not found: ${barcodes}"
        [ sample, fastq, bam, barcodes ]
    }
}

def createPreprocessingChannel() {
    return Channel.fromList(readSamplesheet())
}


//


workflow preprocess {

    take:
        ch  // [sample, fastq_dir | null, bam | null, cell_barcodes]

    main:

        // Stage 1, reads -> MT BAM (skipped for BAM inputs). Hand-over: one library-level MT BAM + barcodes
        fastq = ch.filter { it[1] != null }
        bam = ch.filter { it[2] != null }

        ch_mitobam = bam.map { it -> tuple(it[0], it[2], it[3]) }
        has_fastq = readSamplesheet().any { it[1] != null }

        if (params.assay == "MAESTER" && has_fastq) {

            // STAR-solo, MT reads; the index is built only if some sample starts from FASTQ
            prep_resources()
            get_maester_bam(
                fastq.map { it -> tuple(it[0], it[1]) },
                fastq.map { it -> tuple(it[0], it[3]) },
                prep_resources.out.whitelist,
                prep_resources.out.star_index
            )
            ch_mitobam = ch_mitobam.mix(get_maester_bam.out.mitobam)
            genome = prep_resources.out.genome

        } else if (params.assay == "mtscATAC" && has_fastq) {

            // bwa-mem2 on the NUMT-masked genome, MT reads of the listed cells
            prep_genome()
            get_atac_bam(
                fastq.map { it -> tuple(it[0], it[1]) },
                fastq.map { it -> tuple(it[0], it[3]) },
                prep_genome.out.genome
            )
            ch_mitobam = ch_mitobam.mix(get_atac_bam.out.mitobam)
            genome = prep_genome.out.genome

        } else if (has_fastq) {

            error "FASTQ input is not implemented yet for assay ${params.assay}: give a MT BAM"

        } else {

            prep_genome()
            genome = prep_genome.out.genome

        }

        // Stage 2, MT BAM -> single-cell tables (tool-specific)
        process_mitobam(ch_mitobam, genome)

    emit:

        tables = process_mitobam.out.tables

}


//
