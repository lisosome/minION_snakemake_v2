import pandas as pd
import os
import glob 
import re

def get_mito_samples(wildcards):
    if wildcards.sample in mito_samples:
        return os.path.join(BASE_OUT, wildcards.sample, "1.UNALIGNED_DATA", f"{wildcards.sample}_unaligned.bam")
    else:
        pass


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