// build_tree

// Include here
nextflow.enable.dsl = 2
include { CASSIOPEIA } from "./modules/cassiopeia.nf"
include { IQTREE } from "./modules/iqtree.nf"
include { MPBOOT } from "./modules/mpboot.nf"
include { BOOSTER } from "./modules/booster.nf"

// 
 
//----------------------------------------------------------------------------//
// build_tree subworkflow
//----------------------------------------------------------------------------//

//

workflow build_tree {
    
    take: 
        ch_input   

    main: 

        if (params.tree_algorithm == "cassiopeia") {
            CASSIOPEIA(ch_input)
            trees = CASSIOPEIA.out.tree.groupTuple(by: [0,1])
        } else if (params.tree_algorithm == "mpboot") {
            MPBOOT(ch_input)
            trees = MPBOOT.out.tree.groupTuple(by: [0,1])
        } else if (params.tree_algorithm == "iqtree") {
            IQTREE(ch_input)
            trees = IQTREE.out.tree.groupTuple(by: [0,1])
        } else {
            println('Provide valid tracing system option! (e.g., cassiopeia, mpboot, iqtree)')
        }
        BOOSTER(trees)

    emit:

        final_tree = BOOSTER.out.final_tree
 
}
 