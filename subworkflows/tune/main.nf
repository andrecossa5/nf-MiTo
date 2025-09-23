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

        // Validate that required tuning parameters are arrays
        if (!(params.min_n_positive instanceof List)) {
            error "TUNE workflow requires --min_n_positive to be an array (e.g., '[5,10,15]')"
        }
        if (!(params.af_confident_detection instanceof List)) {
            error "TUNE workflow requires --af_confident_detection to be an array (e.g., '[0.01,0.02,0.05]')"
        }
        if (!(params.min_n_confidently_detected instanceof List)) {
            error "TUNE workflow requires --min_n_confidently_detected to be an array (e.g., '[2,5,10]')"
        }
        if (!(params.min_mean_AD_in_positives instanceof List)) {
            error "TUNE workflow requires --min_mean_AD_in_positives to be an array (e.g., '[1.0,1.25,1.5]')"
        }
        if (!(params.min_AD instanceof List)) {
            error "TUNE workflow requires --min_AD to be an array (e.g., '[1,2,3]')"
        }
        if (!(params.bin_method instanceof List)) {
            error "TUNE workflow requires --bin_method to be an array (e.g., '[\"MiTo\",\"vanilla\"]')"
        }
        if (!(params.t_prob instanceof List)) {
            error "TUNE workflow requires --t_prob to be an array (e.g., '[0.5,0.7,0.9]')"
        }
        if (!(params.min_cell_prevalence instanceof List)) {
            error "TUNE workflow requires --min_cell_prevalence to be an array (e.g., '[0.01,0.05,0.1]')"
        } 

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
