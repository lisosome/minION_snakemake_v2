def get_path(wildcards):
    return WG_df.loc[wildcards.sample, "RAW_DATA"]