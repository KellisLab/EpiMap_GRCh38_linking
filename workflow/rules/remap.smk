### https://stackoverflow.com/questions/61390277/snakemake-how-to-make-a-list-of-input-files-based-on-a-previous-rule-that-prod
wildcard_constraints:
    state="[^\s_]+",
    bssid="BSS\d+"

def list_h5ad_per_mark(wildcards):
    with open(checkpoints.list_marks.get().output[0],'r') as f:
        marks = [ct for ct in f.read().split('\n') if len(ct) > 0]
    return expand("resources/{mark}.h5ad", mark=marks)


rule build_chmm_h5ad:
    input: list_h5ad_per_mark
    output: config["local_links"]["chmm_h5ad"]
    params: chromhmm = config["hyper_links"]["chromhmm"],
            states = config["chromhmm_states"]
    log: "log/build_chmm_h5ad.log"
    conda: "../envs/epimap_hg38.yml"
    shell: "workflow/scripts/concat_mark_h5ad.py -i {input} -c {params.chromhmm} -o {output} -s {params.states}"

rule bigWigToPeaks:
    input:
        annot_table = config["local_links"]["annot_table"],
        ccres = config["local_links"]["ccres_hg19"]
    params:
        mark = lambda wildcards: wildcards.mark,
        bssid = lambda wildcards: wildcards.bssid
    output: "resources/{mark}/{bssid}.tab.gz"
    log: "log/bigwig_to_peaks_{mark}_{bssid}.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "workflow/scripts/bigwig_over_ccres.sh {params.bssid} {params.mark} {input.annot_table} {input.ccres} {output}"

def list_bssid_peaks_per_mark(wildcards):
    with open(checkpoints.download_sample_list.get().output[0], 'r') as f:
        sample_list = [ct for ct in f.read().split('\n') if len(ct) > 0]
    return expand("resources/{mark}/{bssid}.tab.gz", bssid=sample_list, mark=wildcards.mark)

def list_links_per_bssid(wildcards):
    # with open(checkpoints.get_chromhmm_states.get().output[0], 'r') as f:
    #     chmm = [x.strip() for x in f.read().split('\n') if len(ct) > 0]
    return expand("resources/link_{bssid}_{state}.bed.gz", state=config["chromhmm_states"], bssid=wildcards.bssid)

def list_bssid_link_files(wildcards):
    with open(checkpoints.download_sample_list.get().output[0], 'r') as f:
        sample_list = [ct for ct in f.read().split('\n') if len(ct) > 0]
    return expand("resources/link_{bssid}.bed.gz", bssid=sample_list)

rule linking_combine:
    input:
        chromhmm = config["local_links"]["chmm_h5ad"],
        links = list_links_per_bssid
    params:
        bssid = lambda wildcards: wildcards.bssid
    output: "resources/link_{bssid}.bed.gz"
    conda: "../envs/epimap_hg38.yml"
    log: "log/linking_combine_{bssid}.log"
    shell: "python3 workflow/scripts/combine_linking.py --bssid {bssid} --chromhmm {input.chromhmm} --linking {links} -o {output}"

rule linking_perchmm:
    input:
          chromhmm = config["local_links"]["chmm_h5ad"],
          gene_loc = config["local_links"]["gene_bed"],
          rna = config["local_links"]["rna"],
#          marks = ["resources/DNase-seq.h5ad","resources/H3K27ac.h5ad", "resources/H3K4me1.h5ad", "resources/H3K4me2.h5ad", "resources/H3K4me3.h5ad", "resources/H3K9ac.h5ad"]
    params: state = lambda wildcards: wildcards.state,
            bssid = lambda wildcards: wildcards.bssid,
            enh_range = config["params"]["linking_range"],
            overlap = config["params"]["enhancer_overlap"],
            cutoff = config["params"]["linking_cutoff"],
            other = config["params"]["linking_params"]
    output: "resources/link_{bssid}_{state}.bed.gz"
    conda: "../envs/epimap_hg38.yml"
    log: "log/linking_perchromhmm_{bssid}_{state}.log"
    shell: "$(type -P time) python3 workflow/scripts/linking.py --bssid {params.bssid} --state {params.state} --chromhmm {input.chromhmm}  --mark resources/DNase-seq.h5ad resources/H3K27ac.h5ad resources/H3K4me1.h5ad resources/H3K4me2.h5ad resources/H3K4me3.h5ad resources/H3K9ac.h5ad -o {output} --tss {input.gene_loc} --overlap {params.overlap} --range {params.enh_range} --cutoff {params.cutoff} --rna {input.rna} {params.other}"

rule build_mark_h5ad:
    input:
          ccre_list = config["local_links"]["ccres_hg38"],
          sample_list = config["local_links"]["sample_list"],
          peaks = list_bssid_peaks_per_mark,
    params:
          column = "mean",
          mark = lambda wildcards: wildcards.mark
    output: "resources/{mark}.h5ad"
    log: "log/build_mark_h5ad_{mark}.log"
    conda: "../envs/epimap_hg38.yml"
    script: "../scripts/aggregate_peaks_to_h5ad.py"


rule liftover_ccres:
    input:
        liftover = config["local_links"]["lift_over"],
        ccres = config["local_links"]["ccres_hg38"]
    output: config["local_links"]["ccres_hg19"]
    log: "log/lift_over_ccres.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "zcat {input.ccres} | cut -f 1,2,3,5 | liftOver /dev/stdin {input.liftover} {output} /dev/null"

rule download_liftover:
    params: config["hyper_links"]["lift_over"]
    output: config["local_links"]["lift_over"]
    log: "log/lift_over.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "wget {params} -O {output}"

rule download_ccres_hg38:
    params: config["hyper_links"]["ccres_hg38"]
    output: config["local_links"]["ccres_hg38"]
    log: "log/ccres_hg38.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "wget {params} -O {output}"

# checkpoint get_chromhmm_states:
#     params: config["chromhmm_states"]
#     conda: "../envs/epimap_hg38_r.yml"
#     output: config["local_links"]["chmm_list"]
#     shell: "echo ${params} | tr -s ' ' '\n' > {output}"

checkpoint download_sample_list:
    params: config["hyper_links"]["sample_list"]
    output: config["local_links"]["sample_list"]
    log: "log/sample_list.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "wget {params} -O {output}"

checkpoint list_marks:
    input: config["local_links"]["annot_table"]
    output: config["local_links"]["mark_list"]
    log: "log/list_marks.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "cut -f 2 {input} | sort | uniq > {output}"

rule download_annot_table:
    params: config["hyper_links"]["annot_table"]
    output: config["local_links"]["annot_table"]
    log: "log/annot_table.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "wget {params} -O {output}"

rule download_gff:
    params: config["hyper_links"]["gff"]
    output: config["local_links"]["gff"]
    log: "log/download_gff.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "wget {params} -O {output}"

rule convert_gff:
    input: config["local_links"]["gff"]
    output: config["local_links"]["gene_bed"]
    log: "log/convert_gff.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "Rscript workflow/scripts/gff3_to_bed.R {input} {output}"

rule download_rna:
    params: config["hyper_links"]["rna"]
    output: config["local_links"]["rna"]
    log: "log/download_rna.log"
    conda: "../envs/epimap_hg38_r.yml"
    shell: "wget {params} -O {output}"
