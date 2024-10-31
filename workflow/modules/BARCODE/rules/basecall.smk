rule barcode_linking:
    input:
        expand("{raw_data_folders}", raw_data_folders = to_basecall)
    output:
        temp(directory(expand(os.path.join(BASE_OUT, "linking", "tmp_{index}"), index = list(range(len(to_basecall))))))
    log:
        os.path.join(config["paths"]["logdir"], "barcode_linking.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "barcode_linking.tsv")
    resources:
        runtime=8640,
        mem_mb=25000,
        slurm_partition="THIN",
        cpus_per_task=1,
        slurm_account="burlo"
    run:
        for ind, el in enumerate(input):
            tmpdir = output[ind]
            errfile = log[0]
            os.makedirs(tmpdir)
            try:
                for file in glob.iglob(os.path.join(el, '**', '*.pod5'), recursive=True):
                    os.symlink(file, os.path.join(tmpdir, os.path.basename(file)))
            except Exception as e:
                with open(errfile, 'a') as ef:
                    ef.write(str(e) + '\n')

rule basecalling:
    input:
        os.path.join(BASE_OUT, "linking", "tmp_{index}")
    output:
        bam = os.path.join(BASE_OUT, "basecall", "tmp_{index}", "all_sample_basecall_tmp_{index}.bam")
    log:
        os.path.join(config["paths"]["logdir"], "tmp_{index}_basecalling.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "tmp_{index}_basecalling_benchmark.tsv")
    params:
        dorado = config["paths"]["dorado"],
        dorado_args = f"""-v --min-qscore 8 -x auto --kit-name {config["kit-name"]} --trim all"""
    resources:
        runtime=1080,
        mem_mb=75000,
        slurm_partition="GPU",
        cpus_per_task=15,
        slurm_extra="'--gpus=1'",
        slurm_account="burlo"
    shell:
        """
        {params.dorado} basecaller sup@v5.0.0 {input} {params.dorado_args} > {output.bam} 2>> {log}
        """

rule demultiplexing:
    input:
        rules.basecalling.output.bam
    output:
        os.path.join(BASE_OUT, "demultiplex", "tmp_{index}", "tmp_{index}_demultiplex.token")
    log:
        os.path.join(config["paths"]["logdir"], "tmp_{index}_demultiplexing.log")
    benchmark:
        os.path.join(config["paths"]["benchmark"], "tmp_{index}_demultiplexing_benchmark.tsv")
    params:
        dorado = config["paths"]["dorado"],
        dorado_args = f"""--emit-summary --no-classify --no-trim --sort-bam"""
    resources:
        runtime=1080,
        mem_mb=75000,
        slurm_partition="GPU",
        cpus_per_task=15,
        slurm_extra="'--gpus=1'",
        slurm_account="burlo"
    shell:
        """
        outdir=$(dirname {output})
        {params.dorado} demux {params.dorado_args} -o ${{outdir}} {input} 2> {log[0]} && touch {output}
        """

#rule sample_bam_creation:
#    input:
#        expand(os.path.join(BASE_OUT, "demultiplex", "tmp_{index}", "tmp_{index}_demultiplex.token"), index = list(range(len(to_basecall)))),
#    output:
#        #token = touch(os.path.join(BASE_OUT, "{subj}", "{subj}_bam_creation.token")) 
#        bam = os.path.join(BASE_OUT, "{subj}", "1.UNALIGNED_DATA","{subj}_unaligned.bam"),
#    log:
#        os.path.join(config["paths"]["logdir"], "{subj}_bam_creation.log")
#    benchmark:
#        os.path.join(config["paths"]["benchmark"], "{subj}_bam_creation_benchmark.tsv")
#    resources:
#        runtime=8640,
#        mem_mb=150000,
#        slurm_partition="THIN",
#        cpus_per_task=15,
#        slurm_account="burlo"
#    run:
#        try:
#            if len(list(BAR_df.loc[wildcards.subj, "BARCODE"])) > 1:
#                sample_raw_data = list(BAR_df.loc[wildcards.subj, "RAW_DATA"].unique())
#                for fol in sample_raw_data:
#                    tmp = f"tmp_{to_basecall.index(fol)}"
#                    demultiplex_folder = os.path.join(BASE_OUT, "demultiplex", tmp)
#                    sam = wildcards.subj
#                    barcodes = list(BAR_df.query(f"SAMPLE_ID == '{sam}' and RAW_DATA == '{fol}'").BARCODE)
#                    bam_regex = [os.path.join(demultiplex_folder, f"*{kitname}_barcode0{x}.bam") if int(x) < 10 else os.path.join(demultiplex_folder, f"*{kitname}_barcode{x}.bam") for x in barcodes]
#                    globbed = [glob.glob(x)[0] for x in bam_regex]
#                    pysam.merge("-O", "BAM", "-@", "15", '-c', '-p', "-o", output.bam, *[str(gl) for gl in globbed])
#            elif len(list(BAR_df.loc[wildcards.subj, "BARCODE"])) == 1:
#                sample_raw_data = BAR_df.loc[wildcards.subj, "RAW_DATA"]
#                tmp = f"tmp_{to_basecall.index(sample_raw_data)}"
#                demultiplex_folder = os.path.join(BASE_OUT, "demultiplex", tmp)
#                barcodes = BAR_df.loc[wildcards.subj, "BARCODE"]
#                if int(barcodes) < 10: 
#                    bam_regex = os.path.join(demultiplex_folder, f"*{kitname}_barcode0{barcodes}.bam")
#                else:
#                    bam_regex = os.path.join(demultiplex_folder, f"*{kitname}_barcode{barcodes}.bam")
#                globbed = glob.glob(bam_regex)
#                cmd = f"""rsync -avP {globbed[0]} {output.bam}"""
#                shell(cmd)
#        except Exception as e:
#                with open(log[0], 'a') as ef:
#                    ef.write(str(e) + '\n')

