configfile: "config.yml"
rule all:
    input:
        "analyses/data/genbanks/A0901.gbk"
rule download_data:
    output:
        "analyses/data/genbanks/A0901.gbk"
    params: 
        link = config["datalink"]
    log: "logs/downloads/SciLife.log"
    shadow: "minimal"
    shell:
        "wget -O tarfile {params.link} -o {log} && tar -xvzf tarfile >> {log} 2>> {log}"

rule divide_gbks:
    output:
        expand("gbks/{strain}.gbk", strain = config["strains"])
    input:
        expand("analyses/data/genbanks/{strain}.gbk", strain = config["strains"])
    log: "logs/python/gbks.log"
    script:
        "divide_gbks.py"
