configfile: "config.yml"
rule all:
    input:
#        expand("data/genbanks/{strain}.gbk", strain = config["strains"]),
        expand("gbks/{strain}_1.gbk", strain = config["strains"]),
        expand("data/fasta/{GH}.{ext}", GH = config["GHs"], ext = ["faa", "fna"])

#rule download_data:
#    output:
#        expand("data/genbanks/{strain}.gbk", strain = config["strains"])
#    params: 
#        link = config["datalink"]
#    log: "logs/downloads/SciLife.log"
#    shadow: "minimal"
#    shell:
#        "wget -O tarfile {params.link} -o {log} \
#        && tar -xvzf tarfile >> {log} 2>> {log} \
#        && mv analyses/data/genbanks/* data/genbanks"

rule divide_gbks:
    output:
        expand("gbks/{strain}_1.gbk", strain = config["strains"])
    input:
        expand("/home/marina/Akunkeei_files/gbff/{strain}_genomic.gbff", strain = config["strains"])
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

rule separate_GHs:
    output:
        expand("data/fasta/{GH}/{GH}_{sufix}.{ext}", GH = config["GHs"], ext = ["faa", "fna"], sufix = ["repset", "subset"]),
        expand("data/fasta/GH70/{type}_{sufix}.{ext}", type = ["GS1", "GS2", "BRS", "short", "NCB"], sufix = ["repset", "subset"], ext = ["faa", "fna"]),
        expand("data/fasta/GH32/{type}_repset.{ext}", type = ["S1", "S2a", "S2b", "S3"], ext = ["faa", "fna"]),
        expand("data/fasta/GH32/{type}_subset.{ext}", type = ["S2a", "S2b", "S3"], ext = ["faa", "fna"])
    input:
        expand("data/fasta/{GH}/{GH}.{ext}", GH = config["GHs"], ext = ["faa", "fna"])
    params:
        GS1 = config["GS1"],
        GS2 = config["GS2"],
        BRS = config["BRS"],
        short = config["short"],
        NCB = config["NCB"],
        S1 = config["S1"],
        S2a = config["S2a"],
        S2b = config["S2b"],
        S3 = config["S3"],
        repr = config["representatives"],
        subset = config["subset"],
        GH70s = ["GS1", "GS2", "BRS", "short", "NCB"],
        GH32s = ["S1", "S2a", "S2b", "S3"]
    log: "logs/python/subtypes.log"
    conda: "biopython_env.yml"
    script:
        "separate_genes.py"
