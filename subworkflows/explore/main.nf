// explore

// Include here
nextflow.enable.dsl = 2
include { VIZ_MT_SPACE } from "./modules/mt_space_viz.nf"
  
//

//----------------------------------------------------------------------------//
// explore subworkflow
//----------------------------------------------------------------------------//

// Util to create AFM channels
def createAFMChannelWithCoverage() {

    // Validate required parameters
    if (!params.afm_input) {
        error "Error: --afm_input is required for AFM processing"
    }
    ch = Channel.fromPath(params.afm_input)
        .splitCsv(header: true)
        .map { row -> [ row.job_id, row.sample, row.afm, row.coverage ] }

    return ch
}


//


workflow explore {
    
    take:
        ch_input

    main: 
        VIZ_MT_SPACE(ch_input)

    emit:
        plots = VIZ_MT_SPACE.out.plots
        
} 
