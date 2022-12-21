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

rule get_neighboring_genes:
    output:
        expand("data/fasta/other_genes/a_kunkeei_{gname}_{type}.{ext}", gname = config["neighbors"], type = config["GHs"], ext = ["faa", "fna"])
    input:
        expand("/home/marina/Akunkeei_files/gbff/{strain}_genomic.gbff", strain = config["strains"]) #I added GH39 to the gbffs
    conda: "biopython_env.yml"
    log: "logs/python/neighbors.log"
    params: gene_names = config["neighbors"], gtypes = config["GHs"], GH70_subset = config["GH70_subset"], GH32_subset = config["GH32_subset"]
    script:
        "get_neighboring_seqs.py"

#This pipeline requires Mafft, iqtree, pal2nal and paml. It also requires Python3 with biopython and pandas.
#The conda environment files provide all required packages except for pal2nal (v14).

rule msa_gtf:
    output:
        GH70 = expand("data/fasta/GH70/{type}_repset.mafft.{ext}", type = ["GS1", "BRS", "GS2"], ext = ["faa", "fna"]),
        GH32 = expand("data/fasta/GH32/{type}_repset.mafft.{ext}", type = ["S1", "S2a", "S3"], ext = ["faa", "fna"])
    input:
        GH70 = expand("data/fasta/GH70/{type}_repset.{ext}", type = ["GS1", "BRS", "GS2"], ext = ["faa", "fna"]),
        GH32 = expand("data/fasta/GH32/{type}_repset.{ext}", type = ["S1", "S2a", "S3"], ext = ["faa", "fna"])
    threads: 2
    conda: "environment.yml"
    log: "logs/mafft/repset_aln.log"
    shell:
        """
        for i in {input.GH70};
        do
        output=$(echo $i | cut -d'.' -f 1).mafft.$(echo $i | cut -d'.' -f 2)
        mafft-linsi --thread {threads} $i > $output 2>> {log};
        done
        for j in {input.GH32};
        do
        output=$(echo $j | cut -d'.' -f 1).mafft.$(echo $j | cut -d'.' -f 2)
        mafft-linsi --thread {threads} $j > $output 2>> {log};
        done
        """ 

rule msa_other:
    output:
        expand("data/fasta/other_genes/a_kunkeei_{gene}_{type}.mafft.{ext}", gene = config["neighbors"], type = config["GHs"], ext = ["faa", "fna"])
    input:
        expand("data/fasta/other_genes/a_kunkeei_{gene}_{type}.{ext}", gene = config["neighbors"], type = config["GHs"], ext = ["faa", "fna"])
    threads: 2
    conda: "environment.yml"
    log: "logs/mafft/neighbors_aln.log"
    shell:
        """
        for i in {input};
        do
        output=$(echo $i | cut -d'.' -f 1).mafft.$(echo $i | cut -d'.' -f 2)
        mafft-linsi --thread {threads} $i > $output 2>> {log};
        done
        """

rule iqtree_gtf:
    output:
        prot_GH70 = expand("data/fasta/GH70/trees/{type}_repset.mafft.faa.treefile", type = ["GS1", "BRS", "GS2"]),
        gene_GH70 = expand("data/fasta/GH70/trees/{type}_repset.mafft.fna.treefile", type = ["GS1", "BRS", "GS2"]),
        prot_GH32 = expand("data/fasta/GH32/trees/{type}_repset.mafft.faa.treefile", type = ["S2a", "S3"]),
        gene_GH32 = expand("data/fasta/GH32/trees/{type}_repset.mafft.fna.treefile", type = ["S2a", "S3"])
    input:
        prot_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.faa", type = ["GS1", "BRS", "GS2"]),
        gene_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.fna", type = ["GS1", "BRS", "GS2"]),
        prot_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.faa", type = ["S2a", "S3"]),
        gene_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.fna", type = ["S2a", "S3"])
    threads: 12
    conda: "environment.yml"
    log: "logs/iqtree/repset_trees.log"
    shell:
#		'iqtree -s {input.prot} -st AA -m LG+C10+F -bb 1000 -alrt 1000 -v > {output.prot} && iqtree -s {input.gene} -st DNA -m GTR+G4+F -bb 1000 -alrt 1000 -v > {output.gene}'
        """
        mkdir -p data/fasta/GH70/trees && mkdir -p data/fasta/GH32/trees
        for i in {input.prot_GH70};
        do
        iqtree -nt {threads} -s $i -st AA -m LG+G4+F -bb 1000 -bnni >> {log}
        done
        for j in {input.prot_GH32};
        do
        iqtree -nt {threads} -s $j -st AA -m LG+G4+F -bb 1000 -bnni >> {log}
        done
        for k in {input.gene_GH70};
        do
        iqtree -nt {threads} -s $k -st DNA -m GTR+G4+F -bb 1000 -bnni >> {log}
        done
        for l in {input.gene_GH32};
        do
        iqtree -nt {threads} -s $l -st DNA -m GTR+G4+F -bb 1000 -bnni >> {log}
        done
        mv data/fasta/GH70/*.f*a.* data/fasta/GH70/trees
        mv data/fasta/GH32/*.f*a.* data/fasta/GH32/trees
        """ 

rule iqtree_other:
    output:
        expand("data/fasta/other_genes/trees/a_kunkeei_{gene}_{type}.mafft.f{l}a.treefile", gene = config["neighbors"], type = config["GHs"], l = ["a", "n"]),
    input:
        prot = expand("data/fasta/other_genes/a_kunkeei_{gene}_{type}.mafft.faa", gene = config["neighbors"], type = config["GHs"]),
        gene = expand("data/fasta/other_genes/a_kunkeei_{gene}_{type}.mafft.fna", gene = config["neighbors"], type = config["GHs"])
    threads: 4
    conda: "environment.yml"
    log: "logs/iqtree/neighbors_tree.log"
    shell:
        """
        mkdir -p data/fasta/other_genes/trees
        for i in {input.prot};
        do
        iqtree -nt {threads} -s $i -st AA -m LG+G4+F -bb 1000 -bnni >> {log}
        done
        for j in {input.gene};
        do
        iqtree -nt {threads} -s $j -st DNA -m GTR+G4+F -bb 1000 -bnni >> {log}
        done
        mv data/fasta/other_genes/*.f*a.* data/fasta/other_genes/trees
        """    

rule pal2nal:
    output:
        expand("data/codons/{types}_codon.pal2nal", types = ["GS1", "GS2", "BRS", "S2a", "S3"])
    input:
        aln1_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.faa", type = ["GS1", "GS2", "BRS"]),
        aln1_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.faa", type = ["S2a", "S3"]),
	aln2_GH70 = expand("data/fasta/GH70/{type}_repset.fna", type = ["GS1", "GS2", "BRS"]),
        aln2_GH32 = expand("data/fasta/GH32/{type}_repset.fna", type = ["S2a", "S3"])
    log: "logs/pal2nal/repset_codons.log"
    shell:
        """
        pal2nal.v14/pal2nal.pl {input.aln1_GH70[0]} {input.aln2_GH70[0]} -output paml -nogap > {output[0]} 2> {log} &&
        pal2nal.v14/pal2nal.pl {input.aln1_GH70[1]} {input.aln2_GH70[1]} -output paml -nogap > {output[1]} 2>> {log} &&
        pal2nal.v14/pal2nal.pl {input.aln1_GH70[2]} {input.aln2_GH70[2]} -output paml -nogap > {output[2]} 2>> {log} &&
        pal2nal.v14/pal2nal.pl {input.aln1_GH32[0]} {input.aln2_GH32[0]} -output paml -nogap > {output[3]} 2>> {log} &&
        pal2nal.v14/pal2nal.pl {input.aln1_GH32[1]} {input.aln2_GH32[1]} -output paml -nogap > {output[4]} 2>> {log}
        """

rule pal2nal_other:
    output:
        expand("data/codons/a_kunkeei_{CDS}_{type}.pal2nal", CDS = config["neighbors"], type = config["GHs"])
    input:
        aln = expand("data/fasta/other_genes/a_kunkeei_{CDS}_{type}.mafft.faa", CDS = config["neighbors"], type = config["GHs"]),
        gene = expand("data/fasta/other_genes/a_kunkeei_{CDS}_{type}.fna", CDS = config["neighbors"], type = config["GHs"])
    log: "logs/pal2nal/neighbors_codons.log"
    shell:
        """
        for a in {input.aln};
        do
        fna=$(echo $a | cut -d'.' -f 1).fna
        output=data/codons/$(echo $a | cut -d'/' -f 4 | cut -d'.' -f 1).pal2nal
        pal2nal.v14/pal2nal.pl $a $fna -output paml -nogap > $output 2>> {log};
        done
        """

rule codeml_exe:
    output:
        summary = expand("results/{type}/{type}_repset.txt", type = ["GS1", "GS2", "BRS", "S2a", "S3"]),
        dN = expand("results/{type}/2ML.dN", type = ["GS1", "GS2", "BRS", "S2a", "S3"]),
        dS = expand("results/{type}/2ML.dS", type = ["GS1", "GS2", "BRS", "S2a", "S3"])
    input:
        codons = expand("data/codons/{types}_codon.pal2nal", types = ["GS1", "GS2", "BRS", "S2a", "S3"]),
        trees_GH70 = expand("data/fasta/GH70/trees/{type}_repset.mafft.faa.treefile", type = ["GS1", "BRS", "GS2"]),
        trees_GH32 = expand("data/fasta/GH32/trees/{type}_repset.mafft.faa.treefile", type = ["S2a", "S3"])
    params:
        outdir = "/home/marina/GH_project/results"
    log: "logs/python/repset_codeml.log"
    conda: "environment.yml"
    script:
        "codeml_biopython.py"
		
rule run_parser:
    output:
        dNdS = expand("results/{type}/dNdS.tsv", type = ["GS1", "GS2", "BRS", "S2a", "S3"]),
        stats = expand("results/{type}/stats.tsv", type = ["GS1", "GS2", "BRS", "S2a", "S3"])
    input:
        txt = expand("results/{type}/{type}_repset.txt", type = ["GS1", "GS2", "BRS", "S2a", "S3"])
    log: "logs/python/repset_outtabs.log"
    script:
        "parse_codeml.py"
