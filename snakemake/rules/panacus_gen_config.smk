from os.path import join, basename

rule:
    input:
        join(config["graph_dir"], "{sample}.gfa")
    output:
        "./working/{sample}.yaml"
    shell:
       "export GRAPH='{input}'; cat ../../resources/panacus_template.yaml | envsubst > {output}"
