// explore

// Include here
nextflow.enable.dsl = 2
include { VIZ_MT_SPACE } from "./modules/mt_space_viz.nf"
  
//

//----------------------------------------------------------------------------//
// explore subworkflow
//----------------------------------------------------------------------------//


//


workflow explore {
    
    take:
        ch_input

    main: 
        VIZ_MT_SPACE(ch_input)

    emit:
        plots = VIZ_MT_SPACE.out.plots
        
} 
