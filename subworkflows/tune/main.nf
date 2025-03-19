// hyper_tuning

// Include here
nextflow.enable.dsl = 2
include { ONESAMPLE } from "./modules/one_sample.nf"
include { SUMMARISE } from "./modules/summarise.nf"

// 

import java.util.UUID
def generateRandomCode() {
    UUID.randomUUID().toString().replaceAll('-', '').take(10) // Generates a 10-character code
}
        
//

//----------------------------------------------------------------------------//
// hyper_tuning subworkflow
//----------------------------------------------------------------------------//

workflow tune {
    
    take:
        ch_input

    main: 

        ch_input = ch_input.map{ it -> tuple(it[1], it[2]) }
                .combine(Channel.fromList(params.min_n_positive))
                .combine(Channel.fromList(params.af_confident_detection))
                .combine(Channel.fromList(params.min_n_confidently_detected))
                .combine(Channel.fromList(params.min_mean_AD_in_positives))
                .combine(Channel.fromList(params.min_AD))
                .combine(Channel.fromList(params.bin_method))
                .combine(Channel.fromList(params.t_prob))
                .combine(Channel.fromList(params.min_cell_prevalence))
                .map{
                    def (a, b, c, d, e, f, g, h, i, l) = it
                    tuple(a, b, c, d, e, f, g, h, i, l, generateRandomCode())
                }
        ONESAMPLE(ch_input)
        SUMMARISE(ONESAMPLE.out.stats.flatMap{ it -> [it[2], it[3]] }.collect() )

    emit:
        summary = SUMMARISE.out.summary
        
} 
