rule panacus:
    input:
        "./working/{sample}.yaml"
    output:
        "./results/{sample}.html"
    log:
        "./logs/panacus/{sample}.log"
    shell:
        "RUST_LOG=info panacus report {input} > {output} 2> {log}"
