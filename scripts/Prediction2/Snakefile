TRIALS = [
    "2025_AYT_Aurora",
    "24Crk_AY2-3",
    "25_Big6_SVREC_SVREC",
    "AWY1_DVPWA_2024",
    "CornellMaster_2025_McGowan",
    "OHRWW_2025_SPO",
    "STP1_2025_MCG",
    "TCAP_2025_MANKS",
    "YT_Urb_25"
]

VCF_FILES = {
    "2025_AYT_Aurora": "GBS_SDSU_2025.fixed.vcf",
    "24Crk_AY2-3": "GBS_UMN_2020.fixed.vcf",
    "25_Big6_SVREC_SVREC": "ThermoFisher_AgriSeq_4K.fixed.vcf",
    "AWY1_DVPWA_2024": "GBS_WSU_2023.vcf",
    "CornellMaster_2025_McGowan": "GBS_Cornell_2024.vcf",
    "OHRWW_2025_SPO": "Wheat3K.fullfixed.vcf",
    "STP1_2025_MCG": "GBS_TAMU_MCG25.fixed.vcf",
    "TCAP_2025_MANKS": "Allegro_V2.vcf",
    "YT_Urb_25": "GBS_UIUC_2024.fullfixed.vcf",
}

rule all:
    input:
        expand("predictathon_submission/{trial}/submission.csv", trial=TRIALS),
        expand("results/cv0_predictions/{trial}.csv", trial=TRIALS),
        expand("results/cv00_predictions/{trial}.csv", trial=TRIALS),
        "predictathon_submission/ALL_DONE.txt"

# Extract sample names
rule extract_samples:
    output:
        "results/samples/{trial}.txt"
    shell:
        "python -m src.extract_new_samples {wildcards.trial}"

# Filter vcf to sample list
rule filter_vcf:
    input:
        vcf=lambda wc: f"data/predictathon/{wc.trial}/genotypes/{VCF_FILES[wc.trial]}",
        samples="results/samples/{trial}.txt"
    output:
        "data/predictathon/{trial}/genotypes/{trial}_filtered.vcf.gz"
    shell:
        "bcftools view -S {input.samples} -Oz -o {output} {input.vcf}"

# Convert filtered vcf to genotype matrix
rule vcf_to_matrix:
    input:
        "data/predictathon/{trial}/genotypes/{trial}_filtered.vcf.gz"
    output:
        "data/processed/{trial}/geno_matrix.csv"
    shell:
        "python -m src.preprocess_genotypes {wildcards.trial}"

# Build global grm
rule build_global_grm:
    input:
        expand("data/processed/{trial}/geno_matrix.csv", trial=TRIALS)
    output:
        "data/processed/GRM_predictathon.npy",
        "data/processed/GRM_predictathon_lines.txt"
    shell:
        "python -m src.build_global_grm_union"

# Train model
rule train_model:
    input:
        geno="data/processed/{trial}/geno_matrix.csv",
        grm="data/processed/GRM_predictathon.npy",
        pheno="data/processed/unified_training_pheno_mapped.csv"
    output:
        "trained_models/{trial}/final_model.joblib"
    shell:
        "python -m src.train_model {wildcards.trial}"

# Generate cv0 predictions
rule cv0_predict:
    input:
        model="trained_models/{trial}/final_model.joblib",
        grm="data/processed/GRM_predictathon.npy",
        pheno="data/processed/unified_training_pheno_mapped.csv"
    output:
        "results/cv0_predictions/{trial}.csv"
    shell:
        "python -m src.cv0_predict {wildcards.trial}"

# Generate cv00 predictions
rule cv00_predict:
    input:
        model="trained_models/{trial}/final_model.joblib",
        grm="data/processed/GRM_predictathon.npy",
        pheno="data/processed/unified_training_pheno_mapped.csv"
    output:
        "results/cv00_predictions/{trial}.csv"
    shell:
        "python -m src.cv00_predict {wildcards.trial}"

# Build submission
rule build_submission:
    input:
        model="trained_models/{trial}/final_model.joblib",
        grm="data/processed/GRM_predictathon.npy",
        pheno="data/processed/unified_training_pheno_mapped.csv",
        geno="data/processed/{trial}/geno_matrix.csv"
    output:
        "predictathon_submission/{trial}/submission.csv"
    shell:
        "python -m src.build_predictathon_submission {wildcards.trial}"

# Final marker
rule done:
    input:
        expand("predictathon_submission/{trial}/submission.csv", trial=TRIALS),
        expand("results/cv0_predictions/{trial}.csv", trial=TRIALS),
        expand("results/cv00_predictions/{trial}.csv", trial=TRIALS)
    output:
        "predictathon_submission/ALL_DONE.txt"
    shell:
        "echo 'done' > {output}"