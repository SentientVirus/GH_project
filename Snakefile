configfile: "config.yml"
rule all:
    input:
#        expand("data/genbanks/{strain}.gbk", strain = config["strains"]),
#        expand("gbks/{strain}_1.gbk", strain = config["strains"]),
#        expand("data/fasta/{GH}/{GH}.{ext}", GH = config["GHs"], ext = ["faa", "fna"]),
#        "plots/tabfiles/87_A0901_gpr.tab"
HOME = "/home/marina"

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
        expand("gbks/{strain}_1.gbk", strain = config["strains"]),
        expand("gbks/{newstrain}_1.gbk", newstrain = config["extra_strains"])
    input:
        expand("{home}/Akunkeei_files/gbff/{strain}_genomic.gbff", strain = config["strains"], home = HOME),
        expand("{home}/Akunkeei_files/gbff/{newstrain}_genomic.gbff", newstrain = config["extra_strains"], home = HOME)
    params: workdir = HOME + "/GH_project/"
    log: "logs/python/gbks.log"
    script:
        "code/01-divide_gbks.py"

rule retrieve_sequences:
    output:
        expand("data/fasta/{GH}/{GH}.{ext}", GH = config["GHs"], ext = ["faa", "fna"])
    input:
        expand("gbks/{strain}_1.gbk", strain = config["strains"] + config["extra_strains"])
    log: "logs/python/GHs.log"
    conda: "biopython_env.yml"
    script:
        "code/02-retrieve_GH70_domains.py" #"retrieve_GH70s.py"

rule retrieve_full_GH70s:
    output:
        expand("data/fasta/GH70/complete_GH70.{ext}", ext = ["faa", "fna"])
    input:
        expand("gbks/{strain}_1.gbk", strain = config["strains"] + config["extra_strains"])
    log: "logs/python/GH70s.log"
    conda: "biopython_env.yml"
    script:
        "code/03-retrieve_full_GH70s.py"

rule separate_GHs:
    output:
        expand("data/fasta/{GH}/{GH}_{sufix}.{ext}", GH = config["GHs"], ext = ["faa", "fna"], sufix = ["repset", "subset", "all"]),
        expand("data/fasta/GH70/{add}{type}_{sufix}.{ext}", add = ["", "complete_"], type = ["GS1", "GS2", "BRS", "short", "NGB"], sufix = ["repset", "subset", "all"], ext = ["faa", "fna"]),
        expand("data/fasta/GH32/{type}_repset.{ext}", type = ["S1", "S2a", "S2b", "S3"], ext = ["faa", "fna"]),
        expand("data/fasta/GH32/{type}_subset.{ext}", type = ["S1", "S2a", "S2b", "S3"], ext = ["faa", "fna"])
    input:
        expand("data/fasta/{GH}/{GH}.{ext}", GH = config["GHs"], ext = ["faa", "fna"]),
        expand("data/fasta/GH70/complete_GH70.{ext}", ext = ["faa", "fna"])
    params:
        GS1 = config["GS1"],
        GS2 = config["GS2"],
        BRS = config["BRS"],
        short = config["short"],
        NGB = config["NGB"],
        S1 = config["S1"],
        S2a = config["S2a"],
        S2b = config["S2b"],
        S3 = config["S3"],
        repr = config["representatives"],
        subset = config["subset"],
        GH70s = ["GS1", "GS2", "BRS", "short", "NGB"],
        GH32s = ["S1", "S2a", "S2b", "S3"]
    log: "logs/python/subtypes.log"
    conda: "biopython_env.yml"
    script:
        "code/04-separate_genes.py"

rule get_neighboring_genes:
    output:
        expand("data/fasta/other_genes/a_kunkeei_{gname}.{ext}", gname = config["neighbors"], ext = ["faa", "fna"])
    input:
        expand("/home/marina/Akunkeei_files/gbff/modified_gbff/{strain}_genomic.gbff", strain = config["extended_repset"])
    conda: "biopython_env.yml"
    log: "logs/python/neighbors.log"
    params: gene_names = config["neighbors"]
    script:
        "code/02.1-get_neighboring_seqs.py"

rule create_GH70_functional:
    output:
        faa = "data/fasta/GH70/GH70_functional_repset.faa",
        fna = "data/fasta/GH70/GH70_functional_repset.fna",
        faa_outgroup = "data/fasta/GH70/GH70_functional_outgroup_repset.faa"#,
        #fna_outgroup = "data/fasta/GH70_functional_repset_outgroup.fna"
    input:
        faa = expand("data/fasta/GH70/{type}_repset.faa", type = ["GS1", "BRS", "GS2", "NGB"]),
        fna = expand("data/fasta/GH70/{type}_repset.fna", type = ["GS1", "BRS", "GS2", "NGB"]),
        outgroup = "outgroups.fasta" #Replace this to create outgroup gene tree as well
    shell:
        """
        > {output.faa}
        > {output.fna}
        > {output.faa_outgroup}
        cat {input.faa} >> {output.faa}
        cat {input.fna} >> {output.fna}
        cat {output.faa} >> {output.faa_outgroup}
        cat {input.outgroup} >> {output.faa_outgroup}
        """

rule create_GH70_functional_all:
    output:
        faa = "data/fasta/GH70/GH70_functional_all.faa",
        fna = "data/fasta/GH70/GH70_functional_all.fna"
    input:
        faa = expand("data/fasta/GH70/{type}_all.faa", type = ["GS1", "BRS", "GS2", "NGB"]),
        fna = expand("data/fasta/GH70/{type}_all.fna", type = ["GS1", "BRS", "GS2", "NGB"])
    shell:
        """
        > {output.faa}
        > {output.fna}
        cat {input.faa} >> {output.faa}
        cat {input.fna} >> {output.fna}
        """

#This pipeline requires Mafft, iqtree, pal2nal and paml. It also requires Python3 with biopython and pandas.
#The conda environment files provide all required packages except for pal2nal (v14).

rule msa_gtf:
    output:
        GH70 = expand("data/fasta/GH70/{type}_repset.mafft.{ext}", type = ["GH70_functional", "GS1", "BRS", "GS2", "NGB"], ext = ["faa", "fna"]),
        GH32 = expand("data/fasta/GH32/{type}_repset.mafft.{ext}", type = ["GH32", "S1", "S2a", "S2b", "S3"], ext = ["faa", "fna"]),
        GH70_all = expand("data/fasta/GH70/GH70_functional_all.mafft.{ext}", ext = ["faa", "fna"]),
        GH32_all = expand("data/fasta/GH32/GH32_all.mafft.{ext}", ext = ["faa", "fna"]),
        GH70_out = "data/fasta/GH70/GH70_functional_outgroup_repset.mafft.faa"
    input:
        GH70 = expand("data/fasta/GH70/{type}_repset.{ext}", type = ["GH70_functional", "GS1", "BRS", "GS2", "NGB"], ext = ["faa", "fna"]),
        GH32 = expand("data/fasta/GH32/{type}_repset.{ext}", type = ["GH32", "S1", "S2a", "S2b", "S3"], ext = ["faa", "fna"]),
        GH70_all = expand("data/fasta/GH70/GH70_functional_all.{ext}", ext = ["faa", "fna"]),
        GH32_all = expand("data/fasta/GH32/GH32_all.{ext}", ext = ["faa", "fna"]),
        GH70_out = "data/fasta/GH70/GH70_functional_outgroup_repset.faa"
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
        for k in {input.GH70_all};
        do
        output=$(echo $k | cut -d'.' -f 1).mafft.$(echo $k | cut -d'.' -f 2)
        mafft-linsi --thread {threads} $k > $output 2>> {log};
        done
        for l in {input.GH32_all};
        do
        output=$(echo $l | cut -d'.' -f 1).mafft.$(echo $l | cut -d'.' -f 2)
        mafft-linsi --thread {threads} $l > $output 2>> {log};
        done
        mafft-linsi --thread {threads} {input.GH70_out} > {output.GH70_out} 2>> {log};
        """ 

rule msa_other:
    output:
        expand("data/fasta/other_genes/a_kunkeei_{gene}.mafft.{ext}", gene = config["neighbors"], ext = ["faa", "fna"])
    input:
        expand("data/fasta/other_genes/a_kunkeei_{gene}.{ext}", gene = config["neighbors"], ext = ["faa", "fna"])
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
        prot_GH70 = expand("data/fasta/GH70/trees/{type}_repset.mafft.faa.treefile", type = ["GH70_functional", "GS1", "BRS", "GS2", "NGB"]),
        gene_GH70 = expand("data/fasta/GH70/trees/{type}_repset.mafft.fna.treefile", type = ["GH70_functional", "GS1", "BRS", "GS2", "NGB"]),
        prot_GH32 = expand("data/fasta/GH32/trees/{type}_repset.mafft.faa.treefile", type = ["GH32", "S1", "S2a", "S3"]),
        gene_GH32 = expand("data/fasta/GH32/trees/{type}_repset.mafft.fna.treefile", type = ["GH32", "S1", "S2a", "S3"]),
        out_GH70 = "data/fasta/GH70/trees/GH70_functional_outgroup_repset.mafft.faa.treefile"
    input:
        prot_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.faa", type = ["GH70_functional", "GS1", "BRS", "GS2", "NGB"]),
        gene_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.fna", type = ["GH70_functional", "GS1", "BRS", "GS2", "NGB"]),
        prot_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.faa", type = ["GH32", "S1", "S2a", "S3"]),
        gene_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.fna", type = ["GH32", "S1", "S2a", "S3"]),
        out_GH70 = "data/fasta/GH70/GH70_functional_outgroup_repset.mafft.faa"
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
        iqtree -nt {threads} -s {input.out_GH70} -st AA -m LG+G4+F -bb 1000 -bnni >> {log}
        mv data/fasta/GH70/*.f*a.* data/fasta/GH70/trees
        mv data/fasta/GH32/*.f*a.* data/fasta/GH32/trees
        """ 

rule iqtree_other:
    output:
        expand("data/fasta/other_genes/trees/a_kunkeei_{gene}.mafft.f{l}a.treefile", gene = config["neighbors"], l = ["a", "n"]),
    input:
        prot = expand("data/fasta/other_genes/a_kunkeei_{gene}.mafft.faa", gene = config["neighbors"]),
        gene = expand("data/fasta/other_genes/a_kunkeei_{gene}.mafft.fna", gene = config["neighbors"])
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
        expand("data/codons/{types}_codon.pal2nal", types = ["GH70", "GS1", "GS2", "BRS", "GH32", "NGB", "S1", "S2a", "S2b", "S3"])
    input:
        aln1_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.faa", type = ["GH70_functional", "GS1", "GS2", "BRS", "NGB"]),
        aln1_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.faa", type = ["GH32", "S1", "S2a", "S2b", "S3"]),
	aln2_GH70 = expand("data/fasta/GH70/{type}_repset.fna", type = ["GH70_functional", "GS1", "GS2", "BRS", "NGB"]),
        aln2_GH32 = expand("data/fasta/GH32/{type}_repset.fna", type = ["GH32", "S1", "S2a", "S2b", "S3"])
    log: "logs/pal2nal/repset_codons.log"
    shell:
        """
        remove="_functional"
        trim="_repset"
        > {log};
        for input in {input.aln1_GH70};
        do
        output_base="${{input//$trim/}}"
        output_base="${{output_base//$remove/}}"
        fna=$(echo $input | cut -d'.' -f 1).fna
        output=data/codons/$(echo $output_base | cut -d'/' -f 4 | cut -d'.' -f 1)_codon.pal2nal
        pal2nal.v14/pal2nal.pl $input $fna -output paml -nogap > $output 2>> {log};
        done
        for input in {input.aln1_GH32};
        do
        output_base="${{input//$trim/}}"
        fna=$(echo $input | cut -d'.' -f 1).fna
        output=data/codons/$(echo $output_base | cut -d'/' -f 4 | cut -d'.' -f 1)_codon.pal2nal
        pal2nal.v14/pal2nal.pl $input $fna -output paml -nogap > $output 2>> {log};
        done
        """

rule pal2nal_fna:
    output:
        expand("data/codons/{types}_codon.fna", types = ["GH70", "GS1", "GS2", "BRS", "NGB", "GH32", "S1", "S2a", "S2b", "S3"])
    input:
        aln1_GH70 = expand("data/fasta/GH70/{type}_repset.mafft.faa", type = ["GH70_functional", "GS1", "GS2", "BRS", "NGB"]),
        aln1_GH32 = expand("data/fasta/GH32/{type}_repset.mafft.faa", type = ["GH32", "S1", "S2a", "S2b", "S3"]),
        aln2_GH70 = expand("data/fasta/GH70/{type}_repset.fna", type = ["GH70_functional", "GS1", "GS2", "BRS", "NGB"]),
        aln2_GH32 = expand("data/fasta/GH32/{type}_repset.fna", type = ["GH32", "S1", "S2a", "S2b", "S3"])
    log: "logs/pal2nal/repset_fna_codons.log"
    shell:
        """
        remove="_functional"
        trim="_repset"
        > {log};
        for input in {input.aln1_GH70};
        do
        output_base="${{input//$trim/}}"
        output_base="${{output_base//$remove/}}"
        fna=$(echo $input | cut -d'.' -f 1).fna
        output=data/codons/$(echo $output_base | cut -d'/' -f 4 | cut -d'.' -f 1)_codon.fna
        pal2nal.v14/pal2nal.pl $input $fna -output fasta > $output 2>> {log};
        done
        for input in {input.aln1_GH32};
        do
        output_base="${{input//$trim/}}"
        fna=$(echo $input | cut -d'.' -f 1).fna
        output=data/codons/$(echo $output_base | cut -d'/' -f 4 | cut -d'.' -f 1)_codon.fna
        pal2nal.v14/pal2nal.pl $input $fna -output fasta > $output 2>> {log};
        done
        """


rule pal2nal_other:
    output:
        expand("data/codons/a_kunkeei_{CDS}.pal2nal", CDS = config["neighbors"])
    input:
        aln = expand("data/fasta/other_genes/a_kunkeei_{CDS}.mafft.faa", CDS = config["neighbors"]),
        gene = expand("data/fasta/other_genes/a_kunkeei_{CDS}.fna", CDS = config["neighbors"])
    log: "logs/pal2nal/neighbors_codons.log"
    shell:
        """
        > {log};
        for a in {input.aln};
        do
        fna=$(echo $a | cut -d'.' -f 1).fna
        output=data/codons/$(echo $a | cut -d'/' -f 4 | cut -d'.' -f 1).pal2nal
        pal2nal.v14/pal2nal.pl $a $fna -output paml -nogap > $output 2>> {log};
        done
        """

rule codeml_exe:
    output:
        summary = expand("results/{type}/{type}_repset.txt", type = ["GH70", "GS1", "GS2", "BRS", "NGB", "GH32", "S1", "S2a", "S3"]),
        dN = expand("results/{type}/2ML.dN", type = ["GH70", "GS1", "GS2", "BRS", "GH32", "S1", "S2a", "S3"]),
        dS = expand("results/{type}/2ML.dS", type = ["GH70", "GS1", "GS2", "BRS", "GH32", "S1", "S2a", "S3"])
    input:
        codons = expand("data/codons/{types}_codon.pal2nal", types = ["GH70", "GS1", "GS2", "BRS", "NGB", "GH32", "S1", "S2a", "S3"]),
        trees_GH70 = expand("data/fasta/GH70/trees/{type}_repset.mafft.faa.treefile", type = ["GH70_functional", "GS1", "BRS", "GS2", "NGB"]),
        trees_GH32 = expand("data/fasta/GH32/trees/{type}_repset.mafft.faa.treefile", type = ["GH32", "S1", "S2a", "S3"])
    params:
        outdir = "/home/marina/GH_project/results"
    log: "logs/python/repset_codeml.log"
    conda: "environment.yml"
    script:
        "code/11-codeml_biopython.py"

rule codeml_other:
    output:
        summary = expand("results/a_kunkeei_{CDS}/a_kunkeei_{CDS}.txt", CDS = config["neighbors"]),
        dN = expand("results/a_kunkeei_{CDS}/2ML.dN", CDS = config["neighbors"]),
        dS = expand("results/a_kunkeei_{CDS}/2ML.dS", CDS = config["neighbors"])
    input:
        codons = expand("data/codons/a_kunkeei_{CDS}.pal2nal", CDS = config["neighbors"]),
        trees = expand("data/fasta/other_genes/trees/a_kunkeei_{CDS}.mafft.faa.treefile", CDS = config["neighbors"])
    params:
        outdir = "/home/marina/GH_project/results"
    log: "logs/python/neighbors_codeml.log"
    conda: "environment.yml"
    script:
        "code/11.1-codeml_neighbors.py"
		
rule run_parser:
    output:
        dNdS = expand("results/{type}/dNdS.tsv", type = ["GH70", "GS1", "GS2", "BRS", "NGB", "GH32", "S1", "S2a", "S3"]),
        stats = expand("results/{type}/stats.tsv", type = ["GH70", "GS1", "GS2", "BRS", "NGB", "GH32", "S1", "S2a", "S3"])
    input:
        txt = expand("results/{type}/{type}_repset.txt", type = ["GH70", "GS1", "GS2", "BRS", "NGB", "GH32", "S1", "S2a", "S3"])
    log: "logs/python/repset_outtabs.log"
    script:
        "code/12-parse_codeml.py"

rule parse_neighbors:
    output:
        dNdS = expand("results/a_kunkeei_{CDS}/dNdS.tsv", CDS = config["neighbors"]),
        stats = expand("results/a_kunkeei_{CDS}/stats.tsv", CDS = config["neighbors"])
    input:
        txt = expand("results/a_kunkeei_{CDS}/a_kunkeei_{CDS}.txt", CDS = config["neighbors"])
    log: "logs/python/neighbors_outtabs.log"
    script:
        "code/12.1-parse_codeml.py"

## Rules to generate plots
def getTargetFiles():
    targets = []
    for s in config["strains"]:
        no = str(config["no_dict"][s][0])
        while len(no) < 3:
            no = "0" + no
        target = "plots/tabfiles/"+no+"_"+s+"_gpr.tab"
        targets.append(target)
    return targets

output_tabs = getTargetFiles()

rule get_CDS_tabs:
    output:
        output_tabs 
        #lambda wildcards: expand("plots/tabfiles/{no}_{strain}_gpr.tab", strain = config["strains"], no = config["no_dict"][{wildcards.strain}])
    input:
        gbff = expand("gbks/{strain}_1.gbk", strain = config["strains"]),
        tree = "phylogeny.txt"
    log: "logs/python/CDS_tabfiles.log"
#    shadow: "minimal"
    script:
         "code/13-get_CDS_tabs.py"

rule presence_absence_tab:
    output: 
        file = "plots/counts/GH70_32_counts.tab"
    input: 
        tree = "phylogeny.txt",
        GH70 = "data/fasta/GH70/GH70.faa"
    params:
        GS1 = config["GS1"],
        GS2 = config["GS2"],
        BRS = config["BRS"],
        NGB = config["NGB"],
        short = config["short"], 
        S1 = config["S1"],
        S2a = config["S2a"],
        S2b = config["S2b"],
        S3 = config["S3"]
    conda: "biopython_env.yml"
    log: "logs/python/presence_absence.log"
    script:
        "code/14-gtf_CB_per_strain.py"

rule plot_delregion:
    output: "plots/trees/phylogeny.png"
    input: 
        tree = "kunkeei_nonclosed.tree",
        tabs = output_tabs,
        counts = "plots/counts/GH70_32_counts.tab"
    params:
        GS1 = config["GS1"],
        GS2 = config["GS2"],
        BRS = config["BRS"],
        NGB = config["NGB"],
        short = config["short"],
        S1 = config["S1"],
        S2a = config["S2a"],
        S2b = config["S2b"],
        S3 = config["S3"]
    conda: "biopython_env.yml"
    log: "logs/python/plot_delregion.log"
    script: "code/15-ete3_delregion_plot.py"

##Supplementary table
rule suppl_tab:
    output:
        expand("tables/{GH_type}_suppl.tab", GH_type = ["GS1", "GS2", "BRS", "S1", "S2a", "S3", "NGB", "short"])
    params:
        GS1 = config["GS1"],
        GS2 = config["GS2"],
        BRS = config["BRS"],
        NGB = config["NGB"],
        short = config["short"],
        S1 = config["S1"],
        S2a = config["S2a"],
        S2b = config["S2b"],
        S3 = config["S3"],
        representatives = config["representatives"]
    conda: "biopython_env.yml"
    log: "logs/python/suppl_tab.log"
    script: "code/16-get_suppl_tabs.py"

##Add the remaining scripts
