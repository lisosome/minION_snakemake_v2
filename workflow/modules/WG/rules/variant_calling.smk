rule duet:
    input:
        rules.bam_merging.output.bam
    output:
        vcf_sv = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "3.VARIANT.CALLING", "{sample}_sv.vcf.gz"),
        tbi_sv = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "3.VARIANT.CALLING", "{sample}_sv.vcf.gz.tbi"),
        vcf_snp = os.path.join(BASE_OUT,"WG_analyses",  "{sample}", "3.VARIANT.CALLING", "{sample}_phased_snp.vcf.gz"),
        vcf_tbi = os.path.join(BASE_OUT,"WG_analyses",  "{sample}", "3.VARIANT.CALLING", "{sample}_phased_snp.vcf.gz.tbi")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_duet.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_duet.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_duet_benchmark.tsv")
    params:
        ref_genome = config["paths"]["ref_genome"]
    resources:
        runtime=8640,
        mem_mb=100000,
        slurm_partition="EPYC",
        cpus_per_task=50,
        slurm_account="burlo"
    singularity: config["paths"]["wg_container"]
    shell:
        """
        unset CONDA_PREFIX
        export CONDA_PREFIX=/opt2/conda
        
        base=$(dirname {output[0]})
        temp=${{base}}/tmp

        mkdir -p ${{temp}}

        /opt2/conda/bin/duet {input} {params.ref_genome} ${{temp}} -t 50 1> {log[0]} 2> {log[1]} &&
        bcftools view -h ${{temp}}/phased_sv.vcf > ${{temp}}/header.txt 2>> {log[1]} &&
        echo "VALUE {wildcards.sample}" > ${{temp}}/sample_rename.txt 2>> {log[1]} &&
        sed -i '5i ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">' ${{temp}}/header.txt &&
        bgzip ${{temp}}/phased_sv.vcf && tabix -p vcf ${{temp}}/phased_sv.vcf.gz &&
        bcftools reheader -h ${{temp}}/header.txt -s ${{temp}}/sample_rename.txt ${{temp}}/phased_sv.vcf.gz | bcftools sort -Oz -o {output[0]} 1>> {log[0]} 2>> {log[1]} &&
        bcftools index -f -t --threads 50 {output[0]} 1>> {log[0]} 2>> {log[1]} &&
        echo "SAMPLE {wildcards.sample}" > ${{temp}}/sample_rename_snp.txt 2>> {log[1]} &&
        bcftools concat --threads 50 -a ${{temp}}/snp_phasing/phased_chr*.vcf.gz | bcftools sort -Oz -o ${{temp}}/snp_phasing/phased_chrALL.vcf.gz 1>> {log[0]} 2>> {log[1]} && tabix -p vcf ${{temp}}/snp_phasing/phased_chrALL.vcf.gz &&
        bcftools reheader -s ${{temp}}/sample_rename_snp.txt ${{temp}}/snp_phasing/phased_chrALL.vcf.gz | bcftools sort -Oz -o {output[2]} 1>> {log[0]} 2>> {log[1]} &&
        bcftools index -f -t --threads 50 {output[2]} 1>> {log[0]} 2>> {log[1]} &&
        
        rm -r ${{temp}}

        #sed '5i ##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">' ${{temp}}/phased_sv.vcf | bcftools sort -Oz -o {output} 1>> {log[0]} 2>> {log[1]} &&
        
        """

#rule svafotate:
#    input:
#        rules.duet.output.vcf_sv
#    output:
#        vcf = os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "{sample}_svafotate.vcf.gz"),
#        tbi = os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "{sample}_svafotate.vcf.gz.tbi")
#    log:
#        os.path.join(config["paths"]["logdir"], "{sample}_svafotate.log"),
#        os.path.join(config["paths"]["logdir"], "{sample}_svafotate.err")
#    benchmark:
#        os.path.join(config["paths"]["benchmark"], "{sample}_svafotate_benchmark.tsv")
#    resources:
#        runtime=8640,
#        mem_mb=40000,
#        slurm_partition="THIN",
#        cpus_per_task=8,
#        slurm_account="burlo"
#    singularity: "docker://jxprismdocker/prism_svafotate"
#    params:
#        bed = config["paths"]["svafotate_bed"],
#        args = "-f 0.5 --cpu 24 -a best EUR -O vcfgz"
#    shell:
#        """
#        svafotate annotate \
#        -v {input} \
#        -b {params.bed} \
#        {params.args} \
#        -o {output.vcf} && bcftools index -f -t {output.vcf}
#        """

rule all_calls:
    input:
        #rules.svafotate.output.vcf_sv,
        rules.duet.output.vcf_sv,
        rules.duet.output.vcf_snp
    output:
        vcf = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "3.VARIANT.CALLING", "{sample}_all_calls.vcf.gz"),
        tbi = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "3.VARIANT.CALLING", "{sample}_all_calls.vcf.gz.tbi")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_all_calls.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_all_calls.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_all_calls_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    envmodules:
        "bcftools/1.17"
    shell:
        """
        bcftools concat --threads 8 -a {input[0]} {input[1]} | bcftools sort -Oz -o {output[0]} 1> {log[0]} 2> {log[1]} &&
        bcftools index -f -t --threads 8 {output[0]} 1>> {log[0]} 2>> {log[1]}
        """

