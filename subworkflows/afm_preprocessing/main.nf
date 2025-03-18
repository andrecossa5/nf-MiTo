// preprocess

// Include here
nextflow.enable.dsl = 2
include { MITO } from "./modules/mito.nf"
include { scWGS } from "./modules/scWGS.nf"
include { CAS9 } from "./modules/cas9.nf"
include { DISTANCES } from "./modules/distances.nf"
include { DISTANCE_METRICS } from "./modules/distance_metrics.nf"

//
 

// Util to create AFM channels
def createAFMChannel() {

    ch = Channel.fromPath(params.afm_input)
        .splitCsv(header: true)
        .map { row -> [ row.job_id, row.sample, row.afm ] }

    return ch
}


//


// afm preprocessing subworkflow
workflow afm_preprocess {
    
    take:
        ch_jobs 

    main:
     
        // Process input AFMs, based on input scLT_system
        if (params.scLT_system == "MAESTER" || params.scLT_system == "RedeeM") {
            MITO(ch_jobs)   
            ch_afm = MITO.out.afm
        }
        else if (params.scLT_system == "scWGS") {
            scWGS(ch_jobs)   
            ch_afm = scWGS.out.afm
        }
        else if (params.scLT_system == "Cas9") {
            CAS9(ch_jobs)   
            ch_afm = CAS9.out.afm
        }
        else {
            println('Provide valid tracing system option! (e.g., MAESTER, RedeeM, Cas9, scWGS)')
        } 
        
        // Calculate distances
        replicates = Channel.of( 1..(params.n_boot_replicates-1) ).concat(Channel.of( "observed")) 
        DISTANCES(ch_afm.combine(replicates))
        DISTANCE_METRICS(DISTANCES.out.distances.groupTuple(by: [0,1]))

    emit:

        input = DISTANCES.out.distances
        
} 
