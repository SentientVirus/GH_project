configfile: "config.yml"
rule all:
    input:
        expand("data/genbanks/{strain}.gbk", strain = config["strains"]),
        expand("gbks/{strain}_1.gbk", strain = config["strains"]),
        expand("data/fasta/{GH}.{ext}", GH = config["GHs"], ext = ["faa", "fna"])

rule download_data:
    output:
        expand("data/genbanks/{strain}.gbk", strain = config["strains"])
    params: 
        link = config["datalink"]
    log: "logs/downloads/SciLife.log"
    shadow: "minimal"
    shell:
        "wget -O tarfile {params.link} -o {log} \
        && tar -xvzf tarfile >> {log} 2>> {log} \
        && mv analyses/data/genbanks/* data/genbanks"

rule divide_gbks:
    output:
        expand("gbks/{strain}_1.gbk", strain = config["strains"])
    input:
        expand("data/genbanks/{strain}.gbk", strain = config["strains"])
    log: "logs/python/gbks.log"
    script:
        "divide_gbks.py"

rule retrieve_sequences:
    output:
        expand("data/fasta/{GH}/{GH}.{ext}", GH = config["GHs"], ext = ["faa", "fna"])
    input:
        expand("gbks/{strain}_1.gbk", strain = config["strains"])
    log: "logs/python/GHs.log"
    conda: "biopython_env.yml"
    script:
        "retrieve_GH70s.py"
