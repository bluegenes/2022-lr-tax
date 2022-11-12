import os, sys
import pandas as pd
import glob

configfile: "inputs/long-read-study.conf"

out_dir = config.get('output_dir', 'output.protein_gather')
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")

basename = config.get("basename", 'long-read-study')

sample_info = pd.read_csv(config['sample_info'])
SAMPLES = sample_info["name"].to_list()
# rename reads --> raw_reads for clarity
sample_info.rename(columns={'reads': 'raw_reads'}, inplace=True)
# set name as index for easy access
sample_info.set_index('name', inplace=True)

search_databases = config['search_databases'] # must be dictionary


# check params are in the right format, build alpha-ksize combos and param strings
alphabet_info = config['alphabet_info']
alpha_ksize_scaled=[]
nucl_alpha_ksize_scaled=[]
all_param_str=[]
for alpha, info in alphabet_info.items():
    scaled = info["scaled"]
    ksize = info["ksize"]
    if not isinstance(scaled, list):
        scaled = [scaled]
        config["alphabet_info"][alpha]["scaled"] = scaled
    if not isinstance(ksize, list):
        ksize=[ksize]
        config["alphabet_info"][alpha]["ksize"] = ksize
    # build a parameter for the right combinations
    #alpha_ksize += expand(f"{alpha}-k{{ksize}}", ksize = ksize)
    alpha_ksize_scaled += expand(f"{alpha}-k{{ksize}}-sc{{scaled}}", ksize = ksize, scaled=scaled)
    if alpha in ["nucleotide", "dna", "rna"]:
        nucl_alpha_ksize_scaled += expand(f"{alpha}-k{{ksize}}-sc{{scaled}}", ksize = ksize, scaled=scaled)


onstart:
    print("------------------------------")
    print("taxonomic classification with sourmash gather")
    print("------------------------------")


onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")


rule all:
    input:
        ancient(expand(os.path.join(out_dir, f"{basename}.{{read_type}}.queries.zip"), read_type=['raw_reads'])),
        expand(os.path.join(out_dir, '{gather_type}', f"{basename}.{{aks}}.{{ext}}"), aks=alpha_ksize_scaled, gather_type=['gather'], ext=["gather-pathlist.txt", "gather.lineage_summary.tsv"]),
        expand(os.path.join(out_dir, '{gather_type}', '{sample}.{aks}.gather.kreport.txt'), sample=SAMPLES, aks=alpha_ksize_scaled, gather_type=['gather']),
        expand(os.path.join(out_dir, '{gather_type}', '{sample}.{aks}.gather.genbank.kreport.txt'), sample=SAMPLES, aks=nucl_alpha_ksize_scaled, gather_type=['gather']),


rule sourmash_sketch_translate:
    input: ancient(lambda w: sample_info.at[w.sample, w.read_type]) # get raw reads or abundtrim reads. read_type = 'reads' or 'abundtrim'
    output:
        os.path.join(out_dir, "{read_type}", "{sample}.translate.sig.zip")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    log: os.path.join(logs_dir, "sketch_{read_type}", "{sample}.sketch_translate.log")
    benchmark: os.path.join(benchmarks_dir, "sketch_{read_type}", "{sample}.sketch_translate.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch translate {input} -p k=7,k=10,protein,scaled=200,abund \
                                          -p k=16,k=19,dayhoff,scaled=200,abund \
                                          --name {wildcards.sample} -o {output} 2> {log}
        """

rule sourmash_sketch_dna:
    input: ancient(lambda w: sample_info.at[w.sample, w.read_type]) # get raw reads or abundtrim reads. read_type = 'raw_reads' or 'abundtrim'
    output:
        os.path.join(out_dir, "{read_type}", "{sample}.dna.sig.zip")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    log: os.path.join(logs_dir, "sketch_{read_type}", "{sample}.sketch_dna.log")
    benchmark: os.path.join(benchmarks_dir, "sketch_{read_type}", "{sample}.sketch_dna.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch dna {input} -p k=21,k=31,k=51,dna,scaled=1000,abund \
                                    --name {wildcards.sample} -o {output} 2> {log}
        """

localrules: sig_cat
rule sig_cat:
    input: 
        ancient(expand(os.path.join(out_dir, "{{read_type}}", "{sample}.{sketch_type}.sig.zip"), sample=SAMPLES, sketch_type = ['translate', 'dna'])),
    output:
        zipF=os.path.join(out_dir, f"{basename}.{{read_type}}.queries.zip"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    log: os.path.join(logs_dir, "sig_cat_{read_type}", f"{basename}.sigcat.log")
    benchmark: os.path.join(benchmarks_dir, "sig_cat_{read_type}", f"{basename}.sigcat.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sig cat {input} -o {output} 2> {log}
        """

rule gather_sig_from_zipfile:
    input:
        query_zip=ancient(os.path.join(out_dir, f"{basename}.raw_reads.queries.zip")),
        databases = lambda w: search_databases[f"{w.alphabet}-k{w.ksize}"],
    output:
        gather_csv=os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        gather_txt=os.path.join(out_dir, 'gather', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.txt'),
    params:
        #threshold_bp = config.get('sourmash_database_threshold_bp', '0'),
        alpha_cmd = lambda w: "--" + w.alphabet 
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        time=100000,
        partition="bmh",#"bmh",
    log: os.path.join(logs_dir, "gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.log")
    benchmark: os.path.join(benchmarks_dir, "gather", "{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.benchmark")
    #conda: "conf/env/sourmash.yml"
    #conda: "conf/env/sourmash4.5.yml"
    conda: "conf/env/sourmash-dev.yml"
    shell:
        # touch output to let workflow continue in cases where 0 results are found
        """
        echo "DB(s): {input.databases}"
        echo "DB(s): {input.databases}" > {log}

        sourmash sig grep {wildcards.sample} {input.query_zip} {params.alpha_cmd} \
                 --ksize {wildcards.ksize} | sourmash gather - {input.databases} \
                 -o {output.gather_csv} -k {wildcards.ksize} --scaled {wildcards.scaled} \
                 {params.alpha_cmd} > {output.gather_txt} 2>> {log}
        
        touch {output.gather_txt}
        touch {output.gather_csv}
        """
                 #--threshold-bp {params.threshold_bp}  \
                #--threshold-bp={params.threshold_bp}--picklist {input.prefetch_csv}::prefetch \ \ 
                 #--picklist {input.query_picklist}:name:identprefix:exclude \

rule tax_annotate:
    input:
        gather = os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.with-lineages.csv'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'{w.gather_type}'),
        lingather= lambda w: os.path.join(out_dir, f'{w.gather_type}', f'{w.sample}.{w.alphabet}-k{w.ksize}-sc{w.scaled}.gather.with-lineages.csv'),
    #conda: "conf/env/sourmash.yml"
    conda: "conf/env/sourmash4.5.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax annotate -g {input.gather} -t {input.lineages} -o {params.outd}
        """


rule tax_metagenome:
    input:
        gather = os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.krona.tsv'),
        os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.summarized.csv'),
        os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.kreport.txt'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "bml",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'{w.gather_type}'),
        out_base= lambda w: f'{w.sample}.{w.alphabet}-k{w.ksize}-sc{w.scaled}.gather',
    #conda: "conf/env/sourmash.yml"
    #conda: "conf/env/sourmash4.5.yml"
    conda: "conf/env/sourmash-dev.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome -g {input.gather} -t {input.lineages} -o {params.out_base} \
                                --output-dir {params.outd} --output-format krona csv_summary kreport \
                                --rank species
        """


rule tax_metagenome_dna_no_gtdb:
    input:
        gather = os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.csv'),
        lineages = config['database_lineage_files'][:-1], # take off gtdb lineages
    output:
        os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.genbank.krona.tsv'),
        os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.genbank.summarized.csv'),
        os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.genbank.kreport.txt'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "bml",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'{w.gather_type}'),
        out_base= lambda w: f'{w.sample}.{w.alphabet}-k{w.ksize}-sc{w.scaled}.gather.genbank',
    #conda: "conf/env/sourmash.yml"
    #conda: "conf/env/sourmash4.5.yml"
    conda: "conf/env/sourmash-dev.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome -g {input.gather} -t {input.lineages} -o {params.out_base} \
                                --output-dir {params.outd} --output-format krona csv_summary kreport \
                                --rank species
        """


localrules: annotated_gather_csvs_to_pathlist
rule annotated_gather_csvs_to_pathlist:
    input: 
        expand(os.path.join(out_dir, '{{gather_type}}', '{sample}.{{aks}}.gather.with-lineages.csv'), sample=SAMPLES)
    output: 
        os.path.join(out_dir, '{gather_type}', f"{basename}.{{aks}}.gather-pathlist.txt")
    run:
        with open(str(output), 'w') as outF:
            for inF in input:
                outF.write(str(inF)+ "\n")


rule tax_metagenome_lineage_summary:
    input:
        gather_pathlist = os.path.join(out_dir, '{gather_type}', f"{basename}.{{aks}}.gather-pathlist.txt"),
        #gather = expand(os.path.join(out_dir, '{{gather_type}}', '{sample}.{{alphabet}}-k{{ksize}}-sc{{scaled}}.gather.csv'), sample=SAMPLES),
        lineages = config['database_lineage_files'],
    output:
        os.path.join(out_dir, '{gather_type}', f"{basename}.{{aks}}.gather.lineage_summary.tsv")
        #os.path.join(out_dir, '{gather_type}', '{sample}.{alphabet}-k{ksize}-sc{scaled}.gather.lineage_summary.tsv'),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        partition = "low2",
        time=240,
    params:
        outd= lambda w: os.path.join(out_dir, f'{w.gather_type}'),
        out_base= lambda w: f'{basename}.{w.aks}.gather',
    #conda: "conf/env/sourmash.yml"
    conda: "conf/env/sourmash4.5.yml"
    shell:
        """
        mkdir -p {params.outd}
        sourmash tax metagenome --from-file {input.gather_pathlist} -t {input.lineages} -o {params.out_base} --output-dir {params.outd} --output-format lineage_summary --rank species
        """
