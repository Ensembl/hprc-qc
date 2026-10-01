include { PANACUS ; PANACUS_CONFIG } from '../../modules/panacus/main.nf'

workflow PANACUS_STATS {
    take:
    ch_graphs
    config

    main:

    PANACUS_CONFIG(
        ch_graphs,
        file(config.panacus_template),
    )

    PANACUS(
        PANACUS_CONFIG.out,
        ch_graphs,
    )

    emit:
    stats = PANACUS.out
}
