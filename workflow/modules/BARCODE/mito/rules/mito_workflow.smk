rule mt_sample_bam_creation:
    input:
        expand(os.path.join(BASE_OUT, "demultiplex", "tmp_{index}", "tmp_{index}_demultiplex.token"), index = list(range(len(to_basecall)))),
    output:
        #token = touch(os.path.join(BASE_OUT, "{sample}", "{sample}_bam_creation.token")) 
        bam = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "1.UNALIGNED_DATA","{sample}_unaligned.bam"),
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_bam_creation.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_bam_creation_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=150000,
        slurm_partition="THIN",
        cpus_per_task=15,
        slurm_account="burlo"
    run:
        try:
            if len(list(BAR_df.loc[wildcards.sample, "BARCODE"])) > 1:
                sample_raw_data = list(BAR_df.loc[wildcards.sample, "RAW_DATA"].unique())
                for fol in sample_raw_data:
                    tmp = f"tmp_{to_basecall.index(fol)}"
                    demultiplex_folder = os.path.join(BASE_OUT, "demultiplex", tmp)
                    sam = wildcards.sample
                    barcodes = list(BAR_df.query(f"SAMPLE_ID == '{sam}' and RAW_DATA == '{fol}'").BARCODE)
                    bam_regex = [os.path.join(demultiplex_folder, f"*{kitname}_barcode0{x}.bam") if int(x) < 10 else os.path.join(demultiplex_folder, f"*{kitname}_barcode{x}.bam") for x in barcodes]
                    globbed = [glob.glob(x)[0] for x in bam_regex]
                    pysam.merge("-O", "BAM", "-@", "15", '-c', '-p', "-o", output.bam, *[str(gl) for gl in globbed])
            elif len(list(BAR_df.loc[wildcards.sample, "BARCODE"])) == 1:
                sample_raw_data = BAR_df.loc[wildcards.sample, "RAW_DATA"]
                tmp = f"tmp_{to_basecall.index(sample_raw_data)}"
                demultiplex_folder = os.path.join(BASE_OUT, "demultiplex", tmp)
                barcodes = BAR_df.loc[wildcards.sample, "BARCODE"]
                if int(barcodes) < 10: 
                    bam_regex = os.path.join(demultiplex_folder, f"*{kitname}_barcode0{barcodes}.bam")
                else:
                    bam_regex = os.path.join(demultiplex_folder, f"*{kitname}_barcode{barcodes}.bam")
                globbed = glob.glob(bam_regex)
                cmd = f"""rsync -avP {globbed[0]} {output.bam}"""
                shell(cmd)
        except Exception as e:
                with open(log[0], 'a') as ef:
                    ef.write(str(e) + '\n')

rule mito_alignment:
    input:
        #get_mito_samples
        #os.path.join(BASE_OUT, "{sample}", "1.UNALIGNED_DATA","{sample}_unaligned.bam")
        #mito_bam = expand(os.path.join(BASE_OUT, sam, "1.UNALIGNED_DATA", "{sample}_unaligned.bam"), sample=mito_samples),
        rules.mt_sample_bam_creation.output.bam
    output:
        bam = os.path.join(BASE_OUT, "MT_analyses","{sample}", "2.ALIGNMENT","{sample}.bam"),
        bai = os.path.join(BASE_OUT, "MT_analyses","{sample}", "2.ALIGNMENT","{sample}.bam.bai")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_mt_alignment.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_mt_alignment.tsv")
    params:
        dorado = config["paths"]["dorado"],
        mito_ref = config["paths"]["mito_resources"]["ref_mutserve"]
    resources:
        runtime=8640,
        mem_mb=150000,
        slurm_partition="THIN",
        cpus_per_task=15,
        slurm_account="burlo"
    envmodules:
        "samtools/1.17"
    shell:
        """
        {params.dorado} aligner {params.mito_ref} {input} | samtools addreplacerg -r "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" /dev/stdin | samtools sort -O BAM -o {output.bam} 2> {log[0]} &&
        samtools index -@ 15 {output.bam}
        """

rule mito_coverage:
    input:
        rules.mito_alignment.output.bam
    output:
        os.path.join(BASE_OUT, "MT_analyses","{sample}", "2.ALIGNMENT", "coverage", "{sample}.mosdepth.global.dist.txt"),
        os.path.join(BASE_OUT, "MT_analyses","{sample}", "2.ALIGNMENT", "coverage", "{sample}.mosdepth.summary.txt")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_coverage.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_coverage.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_coverage_benchmark.tsv")
    params:
        mosdepth = config["paths"]["mosdepth"],
        mos_opt = "-x -n -t 8",
        out_dir = f"""{BASE_OUT}"""

    resources:
        runtime=8640,
        mem_mb=70000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    shell:
        """
        {params.mosdepth} {params.mos_opt} {params.out_dir}/MT_analyses/{wildcards.sample}/2.ALIGNMENT/coverage/{wildcards.sample} {input} 2> {log[1]} 1> {log[0]}
        """

rule mutserve:
    input:
        rules.mito_alignment.output.bam
    output:
        vcf = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "mutserve", "{sample}.norm.mutserve.vcf.gz"),
        tbi = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "mutserve", "{sample}.norm.mutserve.vcf.gz.tbi")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_mutserve.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_mutserve.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_mutserve_benchmark.tsv")
    singularity: config["paths"]["mito_resources"]["sing_image"]
    params:
        mutserve_opt = "java -jar /opt/mutserve/mutserve.jar call --level 0.01 --mapQ 20 --baseQ 20 --no-ansi --strand-bias 1.6 --write-raw",
        mito_ref = config["paths"]["mito_resources"]["ref_mutserve"],
        baseout = BASE_OUT
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    shell:
        """
        outdir={params.baseout}/MT_analyses/{wildcards.sample}/3.VARIANT.CALLING/mutserve
        mkdir -p ${{outdir}}
        bam_name=$(basename {input})

        echo "${{bam_name}} {wildcards.sample}" > ${{outdir}}/sample_renaming.tsv &&
        {params.mutserve_opt} --reference {params.mito_ref} --output ${{outdir}}/{wildcards.sample}.tmp.mutserve.vcf.gz {input[0]} &&
        bcftools norm \
        -m-any \
        -f {params.mito_ref} \
        -o ${{outdir}}/{wildcards.sample}.tmp.norm.mutserve.vcf.gz -Oz \
        ${{outdir}}/{wildcards.sample}.tmp.mutserve.vcf.gz && 
        bcftools index -t -f ${{outdir}}/{wildcards.sample}.tmp.norm.mutserve.vcf.gz &&
        bcftools reheader -s ${{outdir}}/sample_renaming.tsv ${{outdir}}/{wildcards.sample}.tmp.norm.mutserve.vcf.gz | bcftools view -Oz -o {output.vcf} &&
        bcftools index -f -t {output.vcf}
        """

rule mutect2:
    input:
        rules.mito_alignment.output.bam
    output:
        vcf = os.path.join(BASE_OUT, "MT_analyses","{sample}", "3.VARIANT.CALLING", "mutect2", "{sample}.mutect2.raw.filtered.vcf.gz"),
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_mutect2.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_mutect2.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_mutect2_benchmark.tsv")
    singularity: config["paths"]["mito_resources"]["sing_image"]
    params:
        mutect_opt = "gatk Mutect2 -L chrM --mitochondria-mode --min-base-quality-score 20 -callable-depth 6 --native-pair-hmm-threads 8 --max-reads-per-alignment-start 0",
        mito_ref = config["paths"]["mito_resources"]["ref_mutect"],
        baseout = BASE_OUT
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    shell:
        """
        outdir={params.baseout}/MT_analyses/{wildcards.sample}/3.VARIANT.CALLING/mutect2

        mkdir -p ${{outdir}}

        {params.mutect_opt} -R {params.mito_ref} --tmp-dir ${{outdir}} -I {input[0]} -O {output.vcf}
        """

rule filter_mutect2:
    input:
        rules.mutect2.output.vcf
    output:
        vcf = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "mutect2", "{sample}.norm.mutect2.vcf.gz"),
        tbi = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "mutect2", "{sample}.norm.mutect2.vcf.gz.tbi")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_filter_mutect2.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_filter_mutect2.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_filter_mutect2_benchmark.tsv")
    singularity: config["paths"]["mito_resources"]["sing_image"]
    params:
        mutect_opt = "gatk FilterMutectCalls --min-reads-per-strand 2",
        mito_ref = config["paths"]["mito_resources"]["ref_mutect"],
        baseout = BASE_OUT
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    shell:
        """
        outdir={params.baseout}/MT_analyses/{wildcards.sample}/3.VARIANT.CALLING/mutect2

        {params.mutect_opt} -R {params.mito_ref} --tmp-dir ${{outdir}} -V {input[0]} -O ${{outdir}}/{wildcards.sample}.mutect2.raw.filtered.vcf.gz &&
        bcftools norm -m-any -f {params.mito_ref} -Oz -o {output.vcf} ${{outdir}}/{wildcards.sample}.mutect2.raw.filtered.vcf.gz &&
        bcftools index -f -t {output.vcf}
        """

rule format_mutserve:
    input:
        rules.mutserve.output.vcf
    output:
        filt = os.path.join(BASE_OUT, "MT_analyses","{sample}", "3.VARIANT.CALLING", "filtering", "{sample}.mutserve.filtered.txt")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_filtering_mutserve.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_filtering_mutserve.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_filtering_mutserve_benchmark.tsv")
    params:
        header = "ID\\tFilter\\tPos\\tRef\\tVariant\\tVariantLevel\\tMeanBaseQuality\\tCoverage\\tGT",
        baseout = BASE_OUT
    resources:
        runtime=8640,
        mem_mb=10000,
        slurm_partition="THIN",
        cpus_per_task=2,
        slurm_account="burlo"
    envmodules:
        "bcftools/1.17"
    shell:
        """
        outdir={params.baseout}/MT_analyses/{wildcards.sample}/3.VARIANT.CALLING/filtering
        echo -e "{params.header}" > ${{outdir}}/{wildcards.sample}.mutserve.txt &&
        bcftools query -u -f "{wildcards.sample}\\t%FILTER\\t%POS\\t%REF\\t%ALT\\t[%AF\\t%BQ\\t%DP\\t%GT]\\n" {input[0]} >> ${{outdir}}/{wildcards.sample}.mutserve.txt &&
        awk -F'\\t' 'NR == 1 || (length($4) == 1 && length($5) == 1)' ${{outdir}}/{wildcards.sample}.mutserve.txt > ${{outdir}}/{wildcards.sample}.mutserve.tmp.txt &&
        awk -F'\\t' 'BEGIN {{OFS="\\t"}} {{
            if (NR == 1) {{ print $0, "Type"; next }}
            if ((length($4) > 1 || length($5) > 1) && length($4) != length($5)) {{ $10="3" }}
            else if ($9 == "1") {{ $10="1" }}
            else if ($9 == "0/1" || $9 == "1/0" || $9 == "0|1" || $9 == "1|0") {{ $10="2" }}
            else {{ $10="UNKNOWN" }}
            print
        }}' ${{outdir}}/{wildcards.sample}.mutserve.tmp.txt > {output.filt}
        """

rule format_mutect:
    input:
        rules.mutect2.output.vcf
    output:
        filt = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "filtering", "{sample}.mutect.filtered.txt")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_filtering_mutect.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_filtering_mutect.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_filtering_mutect_benchmark.tsv")
    params:
        header = "ID\\tFilter\\tPos\\tRef\\tVariant\\tVariantLevel\\tMeanBaseQuality\\tCoverage\\tGT",
        baseout = BASE_OUT
    resources:
        runtime=8640,
        mem_mb=10000,
        slurm_partition="THIN",
        cpus_per_task=2,
        slurm_account="burlo"
    envmodules:
        "bcftools/1.17"
    shell:
        """
        outdir={params.baseout}/MT_analyses/{wildcards.sample}/3.VARIANT.CALLING/filtering
        echo -e "{params.header}" > ${{outdir}}/{wildcards.sample}.mutect.txt &&
        bcftools query -u -f "{wildcards.sample}\\t%FILTER\\t%POS\\t%REF\\t%ALT\\t[%AF\\t%BQ\\t%DP\\t%GT]\\n" {input[0]} >> ${{outdir}}/{wildcards.sample}.mutect.txt &&
        awk -F"\\t" 'NR == 1 || ((length($4) > 1 || length($5) > 1) && length($4) != length($5))' ${{outdir}}/{wildcards.sample}.mutect.txt > ${{outdir}}/{wildcards.sample}.mutect.tmp.txt &&
        awk -F'\\t' 'BEGIN {{OFS="\\t"}} {{
            if (NR == 1) {{ print $0, "Type"; next }}
            if ((length($4) > 1 || length($5) > 1) && length($4) != length($5)) {{ $10="3" }}
            else if ($9 == "1") {{ $10="1" }}
            else if ($9 == "0/1" || $9 == "1/0" || $9 == "0|1" || $9 == "1|0") {{ $10="2" }}
            else {{ $10="UNKNOWN" }}
            print
        }}' ${{outdir}}/{wildcards.sample}.mutect.tmp.txt > {output.filt}
        """

rule merging_variants:
    input:
        rules.format_mutserve.output.filt,
        rules.format_mutect.output.filt
    output:
        merged = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "merging", "{sample}.mtDNA.variants.txt")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_merging_variants.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_merging_variants.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_merging_variants_benchmark.tsv")
    singularity: config["paths"]["mito_resources"]["sing_image"]
    params:
        baseout = BASE_OUT,
        sort_opt = "-k ID:N -k Pos:n -k Ref:N -k Type:nr  -k Variant:N"
    resources:
        runtime=8640,
        mem_mb=25000,
        slurm_partition="THIN",
        cpus_per_task=5,
        slurm_account="burlo"
    shell:
        """
        outdir={params.baseout}/MT_analyses/{wildcards.sample}/3.VARIANT.CALLING/merging

        csvtk concat -t {input[0]} {input[1]} --num-cpus 5 -T -o ${{outdir}}/{wildcards.sample}.variants.concat.txt &&
        csvtk sort -t ${{outdir}}/{wildcards.sample}.variants.concat.txt {params.sort_opt} --num-cpus 5 -T -o ${{outdir}}/{wildcards.sample}.variants.sorted.txt &&
        java -jar /opt/VariantMerger.jar ${{outdir}}/{wildcards.sample}.variants.sorted.txt --output {output}
        """

rule annotation:
    input:
        rules.merging_variants.output.merged
    output:
        annotated = os.path.join(BASE_OUT, "MT_analyses","{sample}", "3.VARIANT.CALLING", "annotation", "{sample}.mtDNA.variants.annotated.txt")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_annotation.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_annotation.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_annotation_benchmark.tsv")
    singularity: config["paths"]["mito_resources"]["sing_image"]
    params:
        annot_file = config["paths"]["mito_resources"]["annotation"]
    resources:
        runtime=8640,
        mem_mb=25000,
        slurm_partition="THIN",
        cpus_per_task=5,
        slurm_account="burlo"
    shell:
        """
        java -jar /opt/mutserve/mutserve.jar annotate --input {input[0]} --output {output.annotated} --annotation {params.annot_file}
        """

rule vcf_merging:
    input:
        rules.mutserve.output.vcf,
        rules.filter_mutect2.output.vcf,
        rules.merging_variants.output.merged
    output:
        vcf = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "vcf_merge", "{sample}.filt.SNPs.INDELs.vcf.gz"),
        tbi = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "vcf_merge", "{sample}.filt.SNPs.INDELs.vcf.gz.tbi")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_vcf_merge.log"),
        os.path.join(config["paths"]["logdir"], "{sample}_vcf_merge.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_vcf_merge_benchmark.tsv")
    params:
        baseout = BASE_OUT
    resources:
        runtime=8640,
        mem_mb=10000,
        slurm_partition="THIN",
        cpus_per_task=2,
        slurm_account="burlo"
    envmodules:
        "bcftools/1.17"
    shell:
        """
        outdir={params.baseout}/MT_analyses/{wildcards.sample}/3.VARIANT.CALLING/vcf_merge
        
        bcftools view -v indels {input[1]} -Oz -o ${{outdir}}/{wildcards.sample}.mutect2.indels.vcf.gz && \
        bcftools index -t -f ${{outdir}}/{wildcards.sample}.mutect2.indels.vcf.gz &&

        while read line;do
            echo -e "chrM\\t${{line}}" >> ${{outdir}}/pos_file.tsv
        done < <(tail -n+2 {input[2]} | cut -f3)

        bcftools concat -a {input[0]} ${{outdir}}/{wildcards.sample}.mutect2.indels.vcf.gz | \
        bcftools sort -Oz -o ${{outdir}}/{wildcards.sample}.SNPs_INDELs.vcf.gz && bcftools index -t -f ${{outdir}}/{wildcards.sample}.SNPs_INDELs.vcf.gz &&
        bcftools view -R ${{outdir}}/pos_file.tsv ${{outdir}}/{wildcards.sample}.SNPs_INDELs.vcf.gz -Oz -o {output.vcf} &&
        bcftools index -t -f {output.vcf}
        """

rule mt_vep_annotation:
    input:
        calls=os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "vcf_merge", "{sample}.filt.SNPs.INDELs.vcf.gz"),
        cache=config["vep"]["cache"],
        plugins=config["vep"]["plugin_dir"],
        fasta=config["paths"]["mito_resources"]["ref_mutserve"],
        fai=f"{config["paths"]["mito_resources"]["ref_mutserve"]}.fai",
    output:
        calls=os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "vep", "{sample}.vep_annotated.vcf.gz"),
        stats=os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "vep", "{sample}.variants.html"),
    params:
        plugins=[f"CADD,snv={config.get("vep").get("cadd").get("snv")},indels={config.get("vep").get("cadd").get("indel")}", f"{get_dbNSFP("chrM")}"],
        extra=get_vep_flags("chrM"),
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_vep_annotation.log"),
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_vep_annotation_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=10,
        slurm_account="burlo"
    threads: 10
    wrapper:
        "v4.7.1/bio/vep/annotate"

rule complete_annotation:
    input:
        rules.mt_vep_annotation.output.calls,
        rules.annotation.output.annotated
    output:
        complete = os.path.join(BASE_OUT, "MT_analyses", "{sample}", "3.VARIANT.CALLING", "{sample}.mtDNA.complete_annotation.xlsx")
    log:
        os.path.join(config["paths"]["logdir"], "{sample}_complete_annotation.log"),
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{sample}_complete_annotation_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=25000,
        slurm_partition="THIN",
        cpus_per_task=5,
        slurm_account="burlo"
    run:
        try:
            baseout=os.path.join(BASE_OUT, "MT_analyses", wildcards.sample, "3.VARIANT.CALLING")
            cmd1 = f"""echo "%CHROM,%POS,%Existing_variation,%REF,%ALT,%Consequence,%SYMBOL,%gnomADg_AF,%gnomADg_NFE_AF,%HGVSc,%INTRON,%EXON,%SIFT,%PolyPhen,%CADD_PHRED,%CLIN_SIG,{wildcards.sample}_Genotype" | sed 's/%//g' | sed 's/SYMBOL/GENE/g' | sed 's/CADD_PHRED/CADD/g' | sed 's/CLIN_SIG/ClinVar/g' | sed 's/Existing_variation/rsID/g' > {baseout}/{wildcards.sample}_only_vep.csv"""
            cmd2 = f"""bcftools +split-vep -d -A tab -s worst:any -f "%CHROM,%POS,%Existing_variation,%REF,%ALT,%Consequence,%SYMBOL,%gnomADv3.1_AF,%gnomADv3.1_AF_nfe,%HGVSc,%INTRON,%EXON,%SIFT,%PolyPhen,%CADD_PHRED,%CLIN_SIG[,%GT]\\n" {input[0]} >> {baseout}/{wildcards.sample}_only_vep.csv"""
            shell(f"{cmd1} && {cmd2}")
            
            mut = pd.read_csv(input[1], sep='\t')
            vep = pd.read_csv(f"{baseout}/{wildcards.sample}_only_vep.csv", sep=',')
            col_vec = [40, 41, 42, 43, 44, 45, 39, 6, 7, 8, 9, 11, 13, 14, 20, 21, 22, 23, 24, 28, 25, 26, 27, 30, 31, 32, 33, 47, 48, 52, 53, 54, 55]
            mut['varName'] = mut['Pos'].astype(str) + "_" + mut['Ref'] + "_" + mut['Variant']
            vep['varName'] = vep['POS'].astype(str) + "_" + vep['REF'] + "_" + vep['ALT']
            mut['Type_Label'] = mut['Type'].replace({1: 'Variant', 2: 'Heteroplasmy', 3: 'InDel'})
            mer = pd.merge(mut, vep, on='varName')
            mer = mer.iloc[:, col_vec]
            mer.columns.values[6] = "Variant_type"
            mer.columns.values[7] = "VariantAlleleFraction"
            final = mer.sort_values(by='POS')
            final.to_excel(output.complete, index=False)
        except Exception as e:
            with open(log[0], 'a') as ef:
                ef.write(str(e) + '\n')

#rule mt_cleaning:
#    input:
#        rules.complete_annotation.output.complete
#    output:
#        token = touch(os.path.join(BASE_OUT, "MT_analyses", "{sample}", "{sample}.cleaning.token"))
#    log:
#        os.path.join(config["paths"]["logdir"], "{sample}_cleaning.log"),
#    benchmark:
#        os.path.join(config["paths"]["benchmark"], "{sample}_cleaning_benchmark.tsv")
#    resources:
#        runtime=8640,
#        mem_mb=5000,
#        slurm_partition="THIN",
#        cpus_per_task=1,
#        slurm_account="burlo"
#    params:
#        specific_out = os.path.join(BASE_OUT, "MT_analyses", "{sample}"),
#        baseout = os.path.join(BASE_OUT, "{sample}")
#    shell:
#        """
#        if [[ {params.specific_out}/1.UNALIGNED_DATA/{wildcards.sample}_unaligned.bam ]];then
#            rm -rf {params.baseout}
#        else        
#            mv {params.baseout}/1.UNALIGNED_DATA {params.specific_out} &&
#            rm -rf {params.baseout}
#        fi
#        """