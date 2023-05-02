configfile: "config.yaml"

rule count:
    input:
        config['adapters'],
        # config['twist_sheet'],
        "results/merged/{sb}.fasta"
    