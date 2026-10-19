// Preprocessing subworkflows: reads (FASTQ or BAM) -> single-cell allelic tables
nextflow.enable.dsl = 2

// Subworkflows
include { prep_resources } from "../prep_resources/main"
include { prep_genome } from "../prep_resources/main"
include { get_maester_bam } from "../get_mitobam/main"
include { process_mitobam } from "../process_mitobam/main"


//


def checkParams() {
    // Assay -> pre-processing tool; assays implemented in this version
    def PP_METHODS = ['MAESTER': 'maegatk', 'mtscATAC': 'mgatk', 'ReDeeM': 'redeem-v']
    def IMPLEMENTED = ['MAESTER']
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

        // MAESTER: the STAR index is built only if some sample starts from FASTQ
        ch_mitobam = bam.map { it -> tuple(it[0], it[2], it[3]) }
        if (readSamplesheet().any { it[1] != null }) {
            prep_resources()
            get_maester_bam(
                fastq.map { it -> tuple(it[0], it[1]) },
                fastq.map { it -> tuple(it[0], it[3]) },
                prep_resources.out.whitelist,
                prep_resources.out.star_index
            )
            ch_mitobam = ch_mitobam.mix(get_maester_bam.out.mitobam)
            genome = prep_resources.out.genome
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
