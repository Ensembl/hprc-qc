process PANACUS_CONFIG {
    publishDir 'results/config'

    input:
    path graph
    path base_config

    output:
    path "${graph.getBaseName()}.yaml"

    script:
    """
    cat ${base_config} | sed "s|{{GRAPH}}|${graph}|" > "${graph.getBaseName()}.yaml"
    """
}

process PANACUS {
    conda 'bioconda::panacus'
    publishDir 'results/panacus'

    input:
    path config
    path graph

    output:
    path "${config.getBaseName()}.html"

    script:
    """
	panacus report '${config}' > "${config.getBaseName()}.html"
	"""
}
