import pandas as pd
import os
import glob 
import re

def get_barcodes(to_basecall):
    pairs = {}
    for ind,el in to_basecall:
        kitname = config["kit-name"]
        tmp = f"tmp_{ind}"
        tot_bar = list(BAR_df.query(f"RAW_DATA == {el}").BARCODE)
        file_names = [f"{kitname}_barcode0{x}.bam" if x < 10 else f"{kitname}_barcode0{x}.bam" for x in tot_bar]
        pairs.update({tmp:file_names})


def get_final_output():
    final_output = expand(os.path.join(BASE_OUT, "demultiplex", "tmp_{index}", "tmp_{index}_demultiplex.token"), index = list(range(len(to_basecall))))
    if len(mito_samples) > 0 and len(gs_samples) > 0:
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "1.UNALIGNED_DATA", "{sample}_unaligned.bam"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "2.ALIGNMENT","{sample}.bam"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "2.ALIGNMENT","{sample}.bam.bai"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "2.ALIGNMENT", "coverage", "{sample}.mosdepth.global.dist.txt"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "2.ALIGNMENT", "coverage", "{sample}.mosdepth.summary.txt"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "annotation", "{sample}.mtDNA.variants.annotated.txt"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "vcf_merge", "{sample}.filt.SNPs.INDELs.vcf.gz"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "vcf_merge", "{sample}.filt.SNPs.INDELs.vcf.gz.tbi"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "vep", "{sample}.vep_annotated.vcf.gz"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "vep", "{sample}.variants.html"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "{sample}.mtDNA.complete_annotation.xlsx"), sample = mito_samples))
        #final_output.append("/orfeo/LTS/burlo/LT_storage/shared/resources/vep/indexed_cache")
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "1.UNALIGNED_DATA", "{gsample}_unaligned.bam"), gsample = gs_samples))
        final_output.append(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}""", f"""{config["gene_specific"]["gene_name"]}_for_coverage.bed"""))
        final_output.append(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}""", f"""{config["gene_specific"]["gene_name"]}_ref.fasta"""))
        final_output.append(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}""", f"""{config["gene_specific"]["gene_name"]}_ref.fasta.fai"""))
        final_output.append(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}""", f"""{config["gene_specific"]["gene_name"]}_mask.bed"""))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT","{gsample}.bam"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT","{gsample}.bam.bai"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.global.dist.txt"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.summary.txt"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.region.dist.txt"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.regions.bed.gz"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "3.VARIANT.CALLING", "{gsample}.vcf.gz"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "3.VARIANT.CALLING", "{gsample}.vcf.gz.tbi"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "3.VARIANT.CALLING", "{gsample}.vep_annotated.vcf.gz"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "3.VARIANT.CALLING", "{gsample}.variants.html"), gsample = gs_samples))
    elif len(gs_samples) > 0 and len(mito_samples) == 0:
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "1.UNALIGNED_DATA", "{gsample}_unaligned.bam"), gsample = gs_samples))
        final_output.append(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}""", f"""{config["gene_specific"]["gene_name"]}_for_coverage.bed"""))
        final_output.append(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}""", f"""{config["gene_specific"]["gene_name"]}_ref.fasta"""))
        final_output.append(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}""", f"""{config["gene_specific"]["gene_name"]}_ref.fasta.fai"""))
        final_output.append(os.path.join(BASE_OUT, f"""{config["gene_specific"]["gene_name"]}""", f"""{config["gene_specific"]["gene_name"]}_mask.bed"""))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT","{gsample}.bam"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT","{gsample}.bam.bai"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.global.dist.txt"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.summary.txt"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.region.dist.txt"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "2.ALIGNMENT", "coverage", "{gsample}.mosdepth.regions.bed.gz"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "3.VARIANT.CALLING", "{gsample}.vcf.gz"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "3.VARIANT.CALLING", "{gsample}.vcf.gz.tbi"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "3.VARIANT.CALLING", "{gsample}.vep_annotated.vcf.gz"), gsample = gs_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{gsample}", "3.VARIANT.CALLING", "{gsample}.variants.html"), gsample = gs_samples))
    elif len(mito_samples) > 0 and len(gs_samples) == 0:
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "1.UNALIGNED_DATA", "{sample}_unaligned.bam"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "2.ALIGNMENT","{sample}.bam"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "2.ALIGNMENT","{sample}.bam.bai"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "2.ALIGNMENT", "coverage", "{sample}.mosdepth.global.dist.txt"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "2.ALIGNMENT", "coverage", "{sample}.mosdepth.summary.txt"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "annotation", "{sample}.mtDNA.variants.annotated.txt"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "vcf_merge", "{sample}.filt.SNPs.INDELs.vcf.gz"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "vcf_merge", "{sample}.filt.SNPs.INDELs.vcf.gz.tbi"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "vep", "{sample}.vep_annotated.vcf.gz"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "vep", "{sample}.variants.html"), sample = mito_samples))
        final_output.append(expand(os.path.join(BASE_OUT, "{sample}", "3.VARIANT.CALLING", "{sample}.mtDNA.complete_annotation.xlsx"), sample = mito_samples))
    return final_output


def get_mito_samples(wildcards):
    if wildcards.sample in mito_samples:
        return os.path.join(BASE_OUT, wildcards.sample, "1.UNALIGNED_DATA", f"{wildcards.sample}_unaligned.bam")
    else:
        pass

def get_gs_samples(wildcards):
    if wildcards.gsample in gs_samples:
        return os.path.join(BASE_OUT, wildcards.gsample, "1.UNALIGNED_DATA", f"{wildcards.gsample}_unaligned.bam")
    else:
        pass

def get_mito_bam(wildcards):
    return checkpoints.mito_alignment.get(sample=wildcards.sample).output.bam

def get_gs_bam(wildcards):
    return checkpoints.gs_alignment.get(gsample=wildcards.gsample).output.bam




def get_coord_to_mask(interval):
    chr = interval.split(":")[0]
    gene_start = interval.split(":")[1].split("-")[0]
    gene_end = interval.split(":")[1].split("-")[1]
    return [chr, gene_start, gene_end]


def get_vep_flags(interval):
    extras = ["--af_gnomad", "--af", "--af_1kg", "--variant_class", "--regulatory", "--ccds", "--protein", "--uniprot", "--sift b", "--xref_refseq", "--hgvs", "--hgvsg", "--canonical", "--polyphen b", "--symbol", "--vcf", "--assembly GRCh38"]
    gnomad_annot_string="AC,AN,AF,nhomalt,AC-XY,AN-XY,AF-XY,nhomalt-XY,AC-oth,AN-oth,AF-oth,nhomalt-oth,AC-ami,AN-ami,AF-ami,nhomalt-ami,AC-sas,AN-sas,AF-sas,nhomalt-sas,AC-fin,AN-fin,AF-fin,nhomalt-fin,AC-eas,AN-eas,AF-eas,nhomalt-eas,AC-amr,AN-amr,AF-amr,nhomalt-amr,AC-afr,AN-afr,AF-afr,nhomalt-afr,AC-mid,AN-mid,AF-mid,nhomalt-mid,AC-asj,AN-asj,AF-asj,nhomalt-asj,AC-nfe,AN-nfe,AF-nfe,nhomalt-nfe"
    if "M" in interval:
        gnomad = config.get("vep").get("gnomAD_mito")
        gnomad_str = f"--custom {gnomad},gnomADv3.1,vcf,exact,0,{gnomad_annot_string}"
        extras.append(gnomad_str)
        extra_str = " ".join(extras)
    else:
        gnomad_dir = os.path.dirname(config.get("vep").get("gnomAD_mito"))
        chrom, _, _2 = get_coord_to_mask(interval)
        file_name = f"gnomad.joint.v4.1.sites.{chrom}.vcf.bgz"
        gnomad = os.path.join(gnomad_dir, file_name)
        gnomad_str = f"--custom {gnomad},gnomADv3.1,vcf,exact,0,{gnomad_annot_string}"
        extras.append(gnomad_str)
        extra_str = " ".join(extras)
    return extra_str

def get_dbNSFP(interval):
    dbNSFP_fields="Ensembl_transcriptid,Uniprot_acc,VEP_canonical,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,REVEL_score,REVEL_rankscore,DANN_score,DANN_rankscore,ClinPred_score,ClinPred_rankscore,ClinPred_pred"
    if "M" in interval:
        dbnsfp = os.path.join(config.get("vep").get("dbNSFP_dir"), "dbNSFP4.3a_variant.chrM.gz")
        db_str = f"dbNSFP,{dbnsfp},{dbNSFP_fields}"
    else:
        chrom, _, _2 = get_coord_to_mask(interval)
        file_name = f"dbNSFP4.3a_variant.{chrom}.gz"
        dbnsfp = os.path.join(config.get("vep").get("dbNSFP_dir"), file_name)
        db_str = f"dbNSFP,{dbnsfp},{dbNSFP_fields}"
    return db_str


#def resolve_bams(wildcards):
#    if len(mito_samples) > 0 and len(gs_samples) == 0:
#        return os.path.join(BASE_OUT, wildcards.sample, "1.UNALIGNED_DATA",f"{wildcards.sample}_unaligned.bam")