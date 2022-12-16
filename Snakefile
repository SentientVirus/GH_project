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

rule retrieve_full_GH70s:
    output:
        expand("data/fasta/GH70/complete_GH70.{ext}", ext = ["faa", "fna"])
    input:
        expand("gbks/{strain}_1.gbk", strain = config["strains"])
    log: "logs/python/GH70s.log"
    conda: "biopython_env.yml"
    script:
        "retrieve_full_GH70s.py"

rule separate_GHs:
    output:
        expand("data/fasta/{GH}/{GH}_{sufix}.{ext}", GH = config["GHs"], ext = ["faa", "fna"], sufix = ["repset", "subset"]),
        expand("data/fasta/GH70/{add}{type}_{sufix}.{ext}", add = ["", "complete_"], type = ["GS1", "GS2", "BRS", "short", "NCB"], sufix = ["repset", "subset"], ext = ["faa", "fna"]),
        expand("data/fasta/GH32/{type}_repset.{ext}", type = ["S1", "S2a", "S2b", "S3"], ext = ["faa", "fna"]),
        expand("data/fasta/GH32/{type}_subset.{ext}", type = ["S2a", "S2b", "S3"], ext = ["faa", "fna"])
    input:
        expand("data/fasta/{GH}/{GH}.{ext}", GH = config["GHs"], ext = ["faa", "fna"]),
        expand("data/fasta/GH70/complete_GH70.{ext}", ext = ["faa", "fna"])
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

#This pipeline requires Mafft, iqtree, pal2nal and paml. It also requires Python3 with biopython and pandas.

rule msa_gtf:
    output:
        prot_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.faa", type = ["GS1", "BRS", "GS2"]),
        gene_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.fna", type = ["GS1", "BRS", "GS2"]),
        prot_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.faa", type = ["S1", "S2a", "S3"]),
        gene_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.fna", type = ["S1", "S2a", "S3"])
    input:
        prot_GH70 = expand("data/fasta/GH70/{type}_repset.faa", type = ["GS1", "BRS", "GS2"]),
        gene_GH70 = expand("data/fasta/GH70/{type}_repset.fna", type = ["GS1", "BRS", "GS2"]),
        prot_GH32 = expand("data/fasta/GH32/{type}_repset.faa", type = ["S1", "S2a", "S3"]),
        gene_GH32 = expand("data/fasta/GH32/{type}_repset.fna", type = ["S1", "S2a", "S3"])
    threads: 2
    conda: "environment.yml"
    log: "logs/mafft/repset_aln.log"
    shell:
        """
        mafft-linsi --thread {threads} {input.prot_GH70[0]} > {output.prot_GH70[0]} 2> {log} &&
        mafft-linsi --thread {threads} {input.prot_GH70[1]} > {output.prot_GH70[1]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.prot_GH70[2]} > {output.prot_GH70[2]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.gene_GH70[0]} > {output.gene_GH70[0]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.gene_GH70[1]} > {output.gene_GH70[1]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.gene_GH70[2]} > {output.gene_GH70[2]} 2>> {log} && 
        mafft-linsi --thread {threads} {input.prot_GH32[0]} > {output.prot_GH32[0]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.prot_GH32[1]} > {output.prot_GH32[1]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.prot_GH32[2]} > {output.prot_GH32[2]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.gene_GH32[0]} > {output.gene_GH32[0]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.gene_GH32[1]} > {output.gene_GH32[1]} 2>> {log} &&
        mafft-linsi --thread {threads} {input.gene_GH32[2]} > {output.gene_GH32[2]} 2>> {log}
        """ 

rule iqtree_gtf:
	input:
		prot = 'alignments/{gnrl}_prot.mafft.fasta',
		gene = 'alignments/{gnrl}_genes.mafft.fasta'
	output:
		prot = 'alignments/{gnrl}_prot.mafft.fasta.log',
		gene = 'alignments/{gnrl}_genes.mafft.fasta.log',
		tree = 'alignments/{gnrl}_prot.mafft.fasta.treefile'
	threads: 12
	shell:
#		'iqtree -s {input.prot} -st AA -m LG+C10+F -bb 1000 -alrt 1000 -v > {output.prot} && iqtree -s {input.gene} -st DNA -m GTR+G4+F -bb 1000 -alrt 1000 -v > {output.gene}'
		'iqtree -nt {threads} -s {input.prot} -st AA -m LG+G4+F -bb 1000 -bnni && iqtree -s {input.gene} -st DNA -m GTR+G4+F -bb 1000 -bnni'

rule pal2nal:
	input:
		aln1 = 'alignments/{gnrl}_prot.mafft.fasta',
		aln2 = 'inputs/{gnrl}_genes.fasta'
	output:
		'codons/{gnrl}_codon.pal2nal'
	shell:
		'../pal2nal.v14/pal2nal.pl {input.aln1} {input.aln2} -output paml -nogap > {output}'

rule codeml_exe:
	input:
		'codons/{gnrl}_codon.pal2nal',
		'alignments/{gnrl}_prot.mafft.fasta.treefile'
	output:
		'codeml_out/{gnrl}/{gnrl}.txt',
		'codeml_out/{gnrl}/2ML.dN',
		'codeml_out/{gnrl}/2ML.dS'
	params:
		outdir = '/home/marina/codeml_analysis/subset/codeml_out/{gnrl}'
	script:
		'codeml_biopython.py'
		
rule run_parser:
	input:
		'codeml_out/{gnrl}/{gnrl}.txt'
	output:
		'codeml_out/{gnrl}/dNdS.tsv',
		'codeml_out/{gnrl}/stats.tsv'
	script:
		'../parse_codeml.py'
