def get_reads(sample_df, wildcards, col):
    return sample_df.loc[[wildcards.sample], col].dropna().tolist()


def reads(wildcards):
    if IS_PE:
        return [get_reads(_samples, wildcards, "r1"),
                get_reads(_samples, wildcards, "r2")]
    else:
        return [get_reads(_samples, wildcards, "r1")]


rule preprocess_reads:
    input:
        unpack(reads)
    output:
        os.path.join(config["results"]["preprocess"],
                     "{sample}.interveled.fa")
    log:
        os.path.join(config["logs"]["preprocess"],
                     "{sample}.preprocess.log")
    params:
        prefix = "{sample}"
    run:
        reads_num = len(input)
        if is_PE:
            if reads_num == 2:
                shell("seqtk mergepe %s %s | seqtk seq -A > %s 2> %s" % \
                      (input[0], input[1], output, log))
            else:
                r1_str = " ".join(input[0:reads_num//2-1])
                r2_str = " ".join(input[reads_num//2:])
                r1 = os.path.join(config["results"]["preprocess"],
                                  "%s.merged.1.fq.gz" % params.prefix)
                r1 = os.path.join(config["results"]["preprocess"],
                                  "%s.merged.1.fq.gz" % params.prefix)
                shell("cat %s > %s" % (r1_str, r1))
                shell("cat %s > %s" % (r2_str, r2))
                shell("seqtk mergepe %s %s | seqtk seq -A > %s 2> %s" % \
                      (r1, r2, output, log))
                shell("rm -rf %s %s" % (r1, r2))
        else:
            shell("seqtk seq -A %s > %s 2> %s" % (input[0], output, log))


rule translate_reads:
    input:
        os.path.join(config["results"]["preprocess"],
                     "{sample}.interveled.fa")
    output:
        os.path.join(config["results"]["translate"],
                     "{sample}.pep.fa")
    log:
        os.path.join(config["logs"]["translate"],
                     "{sample}.translate.log")
    params:
        frame = config["params"]["translate"]["frame"],
        codon_table = config["params"]["translate"]["codon_table"]
    shell:
        '''
        transeq {input} {output} \
        -frame {params.frame} \
        -table={params.codon_table} \
        -sformat pearson \
        2> {log}
        '''


def get_sphmms(wildcards):
    return sphmm_df.loc[(wildcards.bgc, wildcards.interval), "sphmm"]


def get_cutoff(wildcards):
    return f1_df.loc[wildcards.bgc, "f1_cutoff"]


# bgc
# sub type
# interval
rule hmm_search:
    input:
        pep = os.path.join(config["results"]["translate"],
                           "{sample}.pep.fa"),
        sphmm = get_sphmm
    output:
        os.path.join(config["results"]["hmm_search"],
                     "{sample}.{bgc}.{interval}.hmm.tbl")
    log:
        os.path.join(config["logs"]["hmm_search"],
                     "{sample}.{bgc}.{interval}.hmm_search.log")
    threads:
        config["params"]["hmm_search"]["threads"]
    params:
        msv_threshold = config["params"]["hmm_search"]["msv_threshold"],
        vit_threshold = config["params"]["hmm_search"]["vit_threshold"],
        fwd_threshold = config["params"]["hmm_search"]["fwd_threshold"]
    shell:
        '''
        hmmsearch --cpu {threads} \
        --F1 {params.msv_threshold} \
        --F2 {params.vit_threshold} \
        --F3 {params.fwd_threshold} \
        --tblout {output} \
        {input.sphmm} {input.pep} > {log}
        '''


rule parse_hmm_results:
    input:
        os.path.join(config["results"]["hmm_search"],
                     "{bgc}/{sample}.{bgc}.{interval}.hmm.tbl")
    output:
        os.path.join(config["results"]["hmm_search"],
                     "{bgc}/{sample}.{bgc}.{interval}.hmm.txt")
    params:
        sample_type = config["params"]["cohort_name"],
        sample_id = "{sample}",
        bgc_type = "{bgc}",
        window = config["params"]["hmm_search"]["window"],
        interval = "{interval}"
    run:
        import pandas as pd
        from metabgc.src.utils import parseHMM, createPandaDF
        result_dict = parseHMM(input, params.sample_type, params.sample_id,
                               params.bgc_type, params.window, params.interval)
       createPandaDF(result_dict, output)


# one bgc, one merge
rule merge_hmm_results:
    input:
        expand(os.path.join(config["results"]["hmm_search"],
                            "{bgc}/{sample}.{bgc}.{{interval}}.hmm.txt"),
               bgc=_bgcs)
    output:
        os.path.join(config["results"]["hmm_search"],
                     "{bgc}/{bgc}.combined.hmm.txt")
    run:
        input_str = " ".join(input)
        shell("cat %s > %s" % (input_str, output))
       

rule identify:
    input:
        hmm = os.path.join(config["results"]["hmm_search"],
                           "{bgc}/{bgc}.combined.hmm.txt"),
        cutoff = get_cutoff
    output:
        filtered = os.path.join(config["results"]["identify"], "{bgc}/spHMM-filtered-results.txt"),
        # os.path.join(config["results"]["bgc_reads"], "{bgc}/{sample}-detected-reads.txt")
        combined = os.path.join(config["results"]["bgc_reads"], "{bgc}/combined.reads_id.txt")
    params:
        identify_dir = os.path.join(config["results"]["identify"], "{bgc}"),
        fasta_dir = os.path.join(config["results"]["bgc_reads"], "{bgc}")
    run:
        import os
        from metabgc.src.metabgcidentify import runidentify
        runidentify(input.hmm, params.identify_dir, input.cutoff, params.fasta_dir)

        shell("cat %s > %s" % (os.path.join(params.fasta_dir, "*-detected-reads.txt"),
                               output.combined))

rule extract:
    input:
        bgc_reads_id = os.path.join(config["results"]["bgc_reads"], "{bgc}/combined.reads_id.txt"),
        samples = expand(os.path.join(config["results"]["preprocess"],
                                      "{sample}.interveled.fa"),
                         sample=_samples.index.unique())
    output:
        os.path.join(config["results"]["bgc_reads"], "{bgc}-identified-biosynthetic-reads.fa")
    run:
        with open(input.bgc_reads_id, 'r') as ih:
            id_list = ih.read()

        with open(output, 'w') as oh:
            for sample_file in input.samples:
                with open(sample_file, 'r') as fh:
                    for seq in SeqIO.parse(fh, 'fasta'):
                        if seq.id in id_list:
                            SeqIO.write(seq, oh, 'fasta')
