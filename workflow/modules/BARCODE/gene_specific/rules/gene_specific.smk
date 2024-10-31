rule gs_sample_bam_creation:
    input:
        expand(os.path.join(BASE_OUT, "demultiplex", "tmp_{index}", "tmp_{index}_demultiplex.token"), index = list(range(len(to_basecall)))),
    output:
        #token = touch(os.path.join(BASE_OUT, "{gsample}", "{gsample}_bam_creation.token")) 
        bam = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "1.UNALIGNED_DATA","{gsample}_unaligned.bam"),
    log:
        os.path.join(config["paths"]["logdir"], "{gsample}_bam_creation.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{gsample}_bam_creation_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=150000,
        slurm_partition="THIN",
        cpus_per_task=15,
        slurm_account="burlo"
    run:
        try:
            if len(list(BAR_df.loc[wildcards.gsample, "BARCODE"])) > 1:
                sample_raw_data = list(BAR_df.loc[wildcards.gsample, "RAW_DATA"].unique())
                for fol in sample_raw_data:
                    tmp = f"tmp_{to_basecall.index(fol)}"
                    demultiplex_folder = os.path.join(BASE_OUT, "demultiplex", tmp)
                    sam = wildcards.gsample
                    barcodes = list(BAR_df.query(f"SAMPLE_ID == '{sam}' and RAW_DATA == '{fol}'").BARCODE)
                    bam_regex = [os.path.join(demultiplex_folder, f"*{kitname}_barcode0{x}.bam") if int(x) < 10 else os.path.join(demultiplex_folder, f"*{kitname}_barcode{x}.bam") for x in barcodes]
                    globbed = [glob.glob(x)[0] for x in bam_regex]
                    pysam.merge("-O", "BAM", "-@", "15", '-c', '-p', "-o", output.bam, *[str(gl) for gl in globbed])
            elif len(list(BAR_df.loc[wildcards.gsample, "BARCODE"])) == 1:
                sample_raw_data = BAR_df.loc[wildcards.gsample, "RAW_DATA"]
                tmp = f"tmp_{to_basecall.index(sample_raw_data)}"
                demultiplex_folder = os.path.join(BASE_OUT, "demultiplex", tmp)
                barcodes = BAR_df.loc[wildcards.gsample, "BARCODE"]
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

rule create_fasta:
    input:
        config["paths"]["ref_genome"]
    output:
        fasta = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""","reference", f"""{config["gene_specific"]["gene_name"]}_ref.fasta"""),
        fai = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""","reference", f"""{config["gene_specific"]["gene_name"]}_ref.fasta.fai"""),
        bed_cov = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""","reference", f"""{config["gene_specific"]["gene_name"]}_for_coverage.bed"""),
        bed_mask = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""","reference", f"""{config["gene_specific"]["gene_name"]}_mask.bed""")
    log:
        os.path.join(config["paths"]["logdir"], "gene_specific_ref.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "gene_specific_ref_benchmark.tsv")
    params:
        chrom = get_coord_to_mask(config["gene_specific"]["coordinates"])[0],
        gene_start = get_coord_to_mask(config["gene_specific"]["coordinates"])[1],
        gene_end = get_coord_to_mask(config["gene_specific"]["coordinates"])[2],
        gene_name = config["gene_specific"]["gene_name"],
        fai = f"{config["paths"]["ref_genome"]}.fai"
    resources:
        runtime=8640,
        mem_mb=150000,
        slurm_partition="THIN",
        cpus_per_task=15,
        slurm_account="burlo"
    envmodules:
        "bedtools2/2.31.1",
        "samtools/1.17"
    shell:
        """
        if echo "{params.chrom}" | grep -q "X"; then
            crom_digit=23
        elif echo "{params.chrom}" | grep -q "Y"; then
            crom_digit=24
        elif echo "{params.chrom}" | grep -q "M"; then
            crom_digit=25
        else
            crom_digit=$(echo "{params.chrom}" | tr -cd '[:digit:]')
        fi
    
        crom_end=$(head -25 {params.fai} | cut -f-2 | awk -F '\\t' 'NR == '"${{crom_digit}}"' {{print $2}}')
        first_end=$(expr {params.gene_start} - 1)
        second_end=$(expr {params.gene_end} + 1)
        outdir=$(dirname {output.fasta})

        echo -e "{params.chrom}\\t1\\t${{first_end}}\\n{params.chrom}\\t${{second_end}}\\t${{crom_end}}" > {output.bed_mask} 2> {log[0]} &&
        echo -e "{params.chrom}\\t{params.gene_start}\\t{params.gene_end}\\t{params.gene_name}" > {output.bed_cov} 2> {log[0]} &&

        bedtools maskfasta -fi <(bedtools getfasta -fi {input} -bed <(echo -e "{params.chrom}\\t1\\t${{crom_end}}")) -bed {output.bed_mask} -fo /dev/stdout | tail -n+2 | cat <(echo ">{params.chrom}") - > {output.fasta} &&
        samtools faidx {output.fasta}
        """

rule gs_alignment:
    input:
        #get_gs_samples,
        rules.gs_sample_bam_creation.output.bam,
        rules.create_fasta.output.fasta,
    output:
        bam = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "2.ALIGNMENT","{gsample}.bam"),
        bai = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "2.ALIGNMENT","{gsample}.bam.bai")
    log:
        os.path.join(config["paths"]["logdir"], "{gsample}_gs_alignment.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{gsample}_gs_alignment.tsv")
    params:
        dorado = config["paths"]["dorado"],
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
        {params.dorado} aligner {input[1]} {input[0]} | samtools addreplacerg -r "@RG\\tID:{wildcards.gsample}\\tSM:{wildcards.gsample}" /dev/stdin | samtools sort -O BAM -o {output.bam} 2> {log[0]} &&
        samtools index -@ 15 {output.bam}
        """

rule gs_coverage:
    input:
        rules.gs_alignment.output.bam,
        rules.create_fasta.output.bed_cov,
    output:
        #os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.global.dist.txt"),
        #os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.summary.txt"),
        #os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.region.dist.txt"),
        #os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.regions.bed.gz"),
        #os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.regions.bed.gz.csi"),
        token = touch(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.token"))
    log:
        os.path.join(config["paths"]["logdir"], "{gsample}_coverage.log"),
        os.path.join(config["paths"]["logdir"], "{gsample}_coverage.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{gsample}_coverage_benchmark.tsv")
    params:
        mosdepth = config["paths"]["mosdepth"],
        mos_opt = "-x -n -t 8",
        out_dir = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""")
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    shell:
        """
        {params.mosdepth} {params.mos_opt} --by {input[1]} {params.out_dir}/{wildcards.gsample}/2.ALIGNMENT/coverage/{wildcards.gsample} {input[0]} 2> {log[1]} 1> {log[0]}
        """

rule gs_clair:
    input:
        rules.gs_alignment.output.bam,
        rules.create_fasta.output.fasta,
        rules.create_fasta.output.bed_cov,
    output:
        vcf = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "merge_output.vcf.gz"),
        tbi = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "merge_output.vcf.gz.tbi"),
        bed = temp(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "variant_calling_bed"))
    log:
        os.path.join(config["paths"]["logdir"], "{gsample}_clair3.log"),
        os.path.join(config["paths"]["logdir"], "{gsample}_clair3.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{gsample}_clair3_benchmark.tsv")
    params:
        clair_opt = """/opt/bin/run_clair3.sh --threads=8 --platform="ont" --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1""",
        clair3_model = config["gene_specific"]["clair3_model"],
        baseout = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""")
    singularity: config["gene_specific"]["clair3_img"]
    envmodules:
        "bcftools/1.17"
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    shell:
        """
        outdir={params.baseout}/{wildcards.gsample}/3.VARIANT.CALLING
        cut -f-3 {input[2]} > {output.bed}
        {params.clair_opt} --bam_fn={input[0]} --ref_fn={input[1]} --bed_fn={output.bed} --model_path={params.clair3_model} --output=${{outdir}}
        """

rule gs_vcf_cleaning:
    input:
        rules.gs_clair.output.vcf,
        rules.create_fasta.output.fasta
    output:
        vcf = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "{gsample}.vcf.gz"),
        tbi = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "{gsample}.vcf.gz.tbi"),
    log:
        os.path.join(config["paths"]["logdir"], "{gsample}_vcf_cleaning.log"),
        os.path.join(config["paths"]["logdir"], "{gsample}_vcf_cleaning.err")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{gsample}_vcf_cleaning_benchmark.tsv")
    envmodules:
        "bcftools/1.17"
    params:
        baseout = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""")
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=8,
        slurm_account="burlo"
    shell:
        """
        outdir={params.baseout}/{wildcards.gsample}/3.VARIANT.CALLING
        bcftools norm -m -any --check-ref -w -f {input[1]} {input[0]} -Oz -o {output.vcf} && bcftools index -f -t {output.vcf} &&
        rm -r ${{outdir}}/log/ ${{outdir}}/tmp/ ${{outdir}}/full* ${{outdir}}/pileup* ${{outdir}}/*.log ${{outdir}}/merge* 
        """

rule gs_vep_annotation:
    input:
        calls=rules.gs_vcf_cleaning.output.vcf,
        cache=config["vep"]["cache"],
        plugins=config["vep"]["plugin_dir"],
        fasta=rules.create_fasta.output.fasta,
        fai=rules.create_fasta.output.fai,
    output:
        calls=os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "{gsample}.vep_annotated.vcf.gz"),
        stats=os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "{gsample}.variants.html"),
    params:
        plugins=[f"CADD,snv={config.get("vep").get("cadd").get("snv")},indels={config.get("vep").get("cadd").get("indel")}", f"{get_dbNSFP(config["gene_specific"]["coordinates"])}"],
        extra=get_vep_flags(config["gene_specific"]["coordinates"]),
    log:
        os.path.join(config["paths"]["logdir"], "{gsample}_vep_annotation.log"),
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{gsample}_vep_annotation_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=40000,
        slurm_partition="THIN",
        cpus_per_task=10,
        slurm_account="burlo"
    threads: 10
    wrapper:
        "v4.7.1/bio/vep/annotate"


rule vep_formatting:
    input:
        rules.gs_vep_annotation.output.calls
    output:
        #csv = temp(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "{gsample}.vep_annotated.csv")),
        excel = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "{gsample}.vep_annotated.xlsx")
    log:
        os.path.join(config["paths"]["logdir"], "{gsample}_vep_formatting.log"),
    benchmark:
        os.path.join(config["paths"]["benchmark"], "{gsample}_vep_formatting_benchmark.tsv")
    resources:
        runtime=8640,
        mem_mb=10000,
        slurm_partition="THIN",
        cpus_per_task=1,
        slurm_account="burlo"
    params:
        csv = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "3.VARIANT.CALLING", "{gsample}.vep_annotated.csv")
    run:
        try:
            cmd1 = f"""echo "%CHROM,%POS,%Existing_variation,%REF,%ALT,%Consequence,%SYMBOL,%gnomADv4.1_AF,%gnomADv4.1_NFE_AF,%HGVSc,%INTRON,%EXON,%SIFT,%PolyPhen,%CADD_PHRED,%CLIN_SIG,{wildcards.gsample}_Genotype" | sed 's/%//g' | sed 's/SYMBOL/GENE/g' | sed 's/CADD_PHRED/CADD/g' | sed 's/CLIN_SIG/ClinVar/g' | sed 's/Existing_variation/rsID/g' > {params.csv}"""
            cmd2 = f"""bcftools +split-vep -d -A tab -s worst:any -f "%CHROM,%POS,%Existing_variation,%REF,%ALT,%Consequence,%SYMBOL,%gnomADv4.1_AF,%gnomADv4.1_AF_nfe,%HGVSc,%INTRON,%EXON,%SIFT,%PolyPhen,%CADD_PHRED,%CLIN_SIG[,%GT]\\n" {input} >> {params.csv}"""
            shell(f"{cmd1} && {cmd2}")
            
            vep = pd.read_csv(params.csv, sep = ",", header=0)
            vep.to_excel(output.excel, index=False)
        except Exception as e:
            with open(log[0], 'a') as ef:
                ef.write(str(e) + '\n')

#rule gs_cleaning:
#    input:
#        rules.vep_formatting.output.excel
#    output:
#        token = touch(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}", "{gsample}.cleaning.token"))
#    log:
#        os.path.join(config["paths"]["logdir"], "{gsample}_cleaning.log"),
#    benchmark:
#        os.path.join(config["paths"]["benchmark"], "{gsample}_cleaning_benchmark.tsv")
#    resources:
#        runtime=8640,
#        mem_mb=5000,
#        slurm_partition="THIN",
#        cpus_per_task=1,
#        slurm_account="burlo"
#    params:
#        specific_out = os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}_analyses""", "{gsample}"),
#        baseout = os.path.join(BASE_OUT, "{gsample}")
#    shell:
#        """
#        if [[ {params.specific_out}/1.UNALIGNED_DATA/{wildcards.gsample}_unaligned.bam ]];then
#            rm -rf {params.baseout}
#        else        
#            mv {params.baseout}/1.UNALIGNED_DATA {params.specific_out} &&
#            rm -rf {params.baseout}
#        fi
#        """