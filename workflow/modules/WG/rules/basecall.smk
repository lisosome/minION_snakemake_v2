rule finding_pod:
    input:
        get_path
    output:
        os.path.join(BASE_OUT,"WG_analyses", "{sample}", "1.BASECALLING", "{sample}_all_pods.tsv")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_all_pods_find.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_all_pods_find_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=25000,
        slurm_partition="THIN",
        cpus_per_task=1,
        slurm_account="burlo"
    shell:
        """
        find {input} -type f -name "*.pod5" > {output} 2> {log}
        """

rule preparing_batches:
    input:
        rules.finding_pod.output
    output:
        scatter.split(os.path.join(BASE_OUT, "WG_analyses", "{{sample}}", "1.BASECALLING", "batch_{scatteritem}", "batch_{scatteritem}.tsv"))
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_batch.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_batch.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_batch_benchmark.tsv")
    #shell:
    #    """
    #    total_lines=$(wc -l {input})
    #    n_lines=$((total_lines/8))
    #    split -l ${{n_lines}} {input} {output} 1> {log[0]} 2> {log[1]}
    #    """
    resources:
        runtime=8640,
        mem_mb=25000,
        slurm_partition="THIN",
        cpus_per_task=1,
        slurm_account="burlo"
    params:
        batches = config["scatter"]["scatter_number"]
    run:        
        with open(input[0], 'r') as f:
            data = f.read().splitlines()

        total_lines = len(data)
        lines_per_file = total_lines // 8
        remainder = total_lines % 8

        start = 0
        for i in range(params.batches):
            end = start + lines_per_file + (1 if i < remainder else 0)
            with open(output[i], 'w') as part_file:
                chunk = data[start:end]
                batch = [f"{x}\n" for x in chunk]
                part_file.writelines(batch)
            start = end

rule linking:
    input:
        os.path.join(BASE_OUT, "WG_analyses", "{sample}","1.BASECALLING", "batch_{scatteritem}", "batch_{scatteritem}.tsv")
    output:
        directory(os.path.join(BASE_OUT, "WG_analyses", "{sample}", "1.BASECALLING", "links_batch_{scatteritem}"))
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_linking_batch_{scatteritem}.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_linking_batch_{scatteritem}.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_linking_batch_{scatteritem}_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=25000,
        slurm_partition="THIN",
        cpus_per_task=1,
        slurm_account="burlo"
    shell:
        """
        mkdir -p {output}
        while read line;do
            ln -s ${{line}} {output}
        done < {input} 1> {log[0]} 2> {log[1]}
        """

rule basecalling:
    input:
        #directory(os.path.join(BASE_OUT, "{{sample}}", "1.BASECALLING", "links_batch_{scatteritem}"))
        rules.linking.output
    output:
        bam = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "2.ALIGNMENT", "batches", "{sample}_{scatteritem}.bam"),
        bai = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "2.ALIGNMENT", "batches", "{sample}_{scatteritem}.bam.bai")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_{scatteritem}_basecalling.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_{scatteritem}_basecalling.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_{scatteritem}_basecalling_benchmark.tsv")
    params:
        dorado = config["paths"]["dorado"],
        dorado_args = f"""-v --min-qscore 8 -x auto --kit-name {config["kit-name"]} --trim all --reference {config["paths"]["ref_genome"]}"""
    resources:
        runtime=1080,
        mem_mb=75000,
        slurm_partition="GPU",
        cpus_per_task=15,
        slurm_extra="'--gpus=1'",
        slurm_account="burlo"
    envmodules:
        "samtools/1.17"
    shell:
        """
        {params.dorado} basecaller sup {input} {params.dorado_args} | samtools addreplacerg -r "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" /dev/stdin | samtools sort -@ 15 -O BAM -o {output[0]} 1> {log[0]} 2> {log[1]} &&
        samtools index {output[0]} 1>> {log[0]} 2>> {log[1]}
        """

rule bam_merging:
    input:
        gather.split(os.path.join(BASE_OUT, "WG_analyses", "{{sample}}", "2.ALIGNMENT", "batches", "{{sample}}_{scatteritem}.bam"))
    output:
        bam = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "2.ALIGNMENT", "{sample}.bam"),
        bai = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "2.ALIGNMENT", "{sample}.bam.bai")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_bam_merging.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_bam_merging.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_bam_merging_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=125000,
        slurm_partition="THIN",
        cpus_per_task=5,
        slurm_account="burlo"
    envmodules:
        "samtools/1.17"
    shell:
        """
        samtools merge -@ 5 -o /dev/stdout {input} | samtools sort -@ 5 -O BAM -o {output[0]} 1> {log[0]} 2> {log[1]} &&
        samtools index -@ 5 {output[0]} 2>> {log[1]}
        """

rule coverage:
    input:
        rules.bam_merging.output.bam
    output:
        os.path.join(BASE_OUT, "WG_analyses", "{sample}", "2.ALIGNMENT", "coverage", "{sample}.mosdepth.global.dist.txt"),
        os.path.join(BASE_OUT, "WG_analyses", "{sample}", "2.ALIGNMENT", "coverage", "{sample}.mosdepth.summary.txt")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_coverage.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_coverage.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_coverage_benchmark.tsv")
    params:
        mosdepth = config["paths"]["mosdepth"],
        mos_opt = "-x -n -t 8",
        out_dir = f"""{BASE_OUT}"""
        #outdir = os.path.join(BASE_OUT, "WG_analyses", "{sample}", "2.ALIGNMENT", "coverage")
    resources:
        runtime=8640,
        mem_mb=70000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    shell:
        """
        {params.mosdepth} {params.mos_opt} {params.out_dir}/WG_analyses/{wildcards.sample}/2.ALIGNMENT/coverage/{wildcards.sample} {input} 2> {log[1]} 1> {log[0]}
        """
