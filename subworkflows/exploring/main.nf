// explore_spaces

// Include here
nextflow.enable.dsl = 2
include { VIZ_MT_SPACE } from "./modules/mt_space_viz.nf"
  
//

//----------------------------------------------------------------------------//
// explore_spaces subworkflow
//----------------------------------------------------------------------------//

workflow explore_spaces {
    
    take:
        ch_input

    main: 
        VIZ_MT_SPACE(ch_input)

    emit:
        plots = VIZ_MT_SPACE.out.plots
        
} 
