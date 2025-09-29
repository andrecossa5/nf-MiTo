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

        // Convert parameters to lists for flexible handling
        def min_n_positive_list = params.min_n_positive instanceof List ? params.min_n_positive : [params.min_n_positive]
        def af_confident_detection_list = params.af_confident_detection instanceof List ? params.af_confident_detection : [params.af_confident_detection]
        def min_n_confidently_detected_list = params.min_n_confidently_detected instanceof List ? params.min_n_confidently_detected : [params.min_n_confidently_detected]
        def min_mean_AD_in_positives_list = params.min_mean_AD_in_positives instanceof List ? params.min_mean_AD_in_positives : [params.min_mean_AD_in_positives]
        def min_AD_list = params.min_AD instanceof List ? params.min_AD : [params.min_AD]
        def bin_method_list = params.bin_method instanceof List ? params.bin_method : [params.bin_method]
        def t_prob_list = params.t_prob instanceof List ? params.t_prob : [params.t_prob]
        def min_cell_prevalence_list = params.min_cell_prevalence instanceof List ? params.min_cell_prevalence : [params.min_cell_prevalence]

        ch_input = ch_input.map{ it -> tuple(it[1], it[2]) }
                .combine(Channel.fromList(min_n_positive_list))
                .combine(Channel.fromList(af_confident_detection_list))
                .combine(Channel.fromList(min_n_confidently_detected_list))
                .combine(Channel.fromList(min_mean_AD_in_positives_list))
                .combine(Channel.fromList(min_AD_list))
                .combine(Channel.fromList(bin_method_list))
                .combine(Channel.fromList(t_prob_list))
                .combine(Channel.fromList(min_cell_prevalence_list))
                .map{
                    def (a, b, c, d, e, f, g, h, i, l) = it
                    tuple(a, b, c, d, e, f, g, h, i, l, generateRandomCode())
                }
        ONESAMPLE(ch_input)
        SUMMARISE(ONESAMPLE.out.stats.flatMap{ it -> [it[2], it[3]] }.collect() )

    emit:
        summary = SUMMARISE.out.summary
        
} 
