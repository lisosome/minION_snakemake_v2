import pandas as pd
import os
import glob 


def get_gs_samples(wildcards):
    if wildcards.gsample in gs_samples:
        return os.path.join(BASE_OUT, wildcards.gsample, "1.UNALIGNED_DATA", f"{wildcards.gsample}_unaligned.bam")
    else:
        pass


def get_coord_to_mask(interval):
    chr = interval.split(":")[0]
    gene_start = interval.split(":")[1].split("-")[0]
    gene_end = interval.split(":")[1].split("-")[1]
    return [chr, gene_start, gene_end]


def get_vep_flags(interval):
    extras = ["--af_gnomad", "--af", "--af_1kg", "--variant_class", "--regulatory", "--ccds", "--protein", "--uniprot", "--sift b", "--xref_refseq", "--hgvs", "--hgvsg", "--canonical", "--polyphen b", "--symbol", "--vcf", "--assembly GRCh38"]
    gnomad_annot_string="AC,AN,AF,nhomalt,AC_XY,AN_XY,AF_XY,nhomalt_XY,AC_oth,AN_oth,AF_oth,nhomalt_oth,AC_ami,AN_ami,AF_ami,nhomalt_ami,AC_sas,AN_sas,AF_sas,nhomalt_sas,AC_fin,AN_fin,AF_fin,nhomalt_fin,AC_eas,AN_eas,AF_eas,nhomalt_eas,AC_amr,AN_amr,AF_amr,nhomalt_amr,AC_afr,AN_afr,AF_afr,nhomalt_afr,AC_mid,AN_mid,AF_mid,nhomalt_mid,AC_asj,AN_asj,AF_asj,nhomalt_asj,AC_nfe,AN_nfe,AF_nfe,nhomalt_nfe"
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
        gnomad_str = f"--custom {gnomad},gnomADv4.1,vcf,exact,0,{gnomad_annot_string}"
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
