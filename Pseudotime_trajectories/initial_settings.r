cols_scale <- colorRampPalette(c("#2C7BB6", "white", "#D7191C"))(100)

cols_class_uni <- c("EN" = "#197EC0",
    "IN" = "#1A9993",
    "Immune" = "#C80813",
    "Mural" = "#FD8CC1",
    "Endo" = "#FD7446",
    "Astro" = "#D2AF81",
    "Oligo" = "#D5E4A2",
    "OPC" = "#FED439",
    "Others" = "#BFBFBF")

cols_subclass_uni <- c("EN_L2_3_IT" = "#659EC7", "EN_L3_5_IT_1" = "#95B9C7", "EN_L3_5_IT_2" = "#6495ED", "EN_L3_5_IT_3" = "#79BAEC", "EN_L5_6_NP" = "#0020C2", "EN_L5_ET" = "#4863A0", "EN_L6_CT" = "#3BB9FF", "EN_L6_IT_1" = "#B4CFEC", "EN_L6_IT_2" = "#1589FF", "EN_L6B" = "#488AC7", "EN_dev" = "#0433FF",
    "MGE_dev" = "#8CCCC9", "CGE_dev" = "#0D4C4A",
    "IN_ADARB2" = "#4E8975", "IN_LAMP5_LHX6" = "#7BCCB5", "IN_LAMP5_RELN" = "#3B9C9C", "IN_PVALB" = "#B2C248", "IN_PVALB_CHC" = "#008080", "IN_SST" = "#728C00", "IN_VIP" = "#89C35C", 
    "Micro" = "#F75D59", "PVM" = "#C11B17", "Adaptive" = "#6F4E37", "PC" = "#E0B0FF", "SMC" = "#4E387E", "VLMC" = "#B93B8F", "Endo" = "#FFA62F", 
    "Astro" = "#C19A6B", 
    "Oligo" = "#ECE5B6", "OPC" = "#FFF380")

cols_traj <- c("Deep-non-IT" = "#0c3f60", "Deep-IT" = "#14659a", "Upper-IT" = "#8cbee0", 
    "L2_3_IT" = "#659EC7", "L3_5_IT_1" = "#95B9C7", "L3_5_IT_2" = "#6495ED", "L3_5_IT_3" = "#79BAEC", "L5_6_NP" = "#0020C2", "L5_ET" = "#4863A0", "L6_CT" = "#3BB9FF", "L6_IT_1" = "#B4CFEC", "L6_IT_2" = "#1589FF", "L6B" = "#488AC7",
    "MGE" = "#8CCCC9", "CGE" = "#0D4C4A", 
    "LAMP5_LHX6" = "#7BCCB5", "LAMP5_RELN" = "#3B9C9C", "PVALB" = "#B2C248", "PVALB_CHC" = "#008080", "SST" = "#728C00", "VIP_TRPC6" = "#89C35C", "VIP_BCL11B" = "#A9D18E",
    "ADARB2_RAB37" = "#58937F", "ADARB2_COL12A1" = "#447F6B", "ADARB2_SYT10" = "#4E9375", "ADARB2_SV2C" = "#4E7F75",
    "Micro" = "#F75D59", 
    "PA" = "#9A7B56", "FA" = "#F4CD9E",
    "Oligo" = "#ECE5B6")

cols_stage_id <- c("Fetal" = "#F4ECD3", "Neonatal" = "#E9D8A6", "Infancy" = "#F1D4A4", "Childhood" = "#F9CFA1", "Adolescence" = "#FFCC66", "Young_Adulthood" = "#EE9B00", "Middle_Adulthood" = "#CA6702", "Late_Adulthood" = "#8C510A")

trait_info <- c(
    "sz3" = "SCZ",
    "bip2" = "BD",
    "mdd_ipsych" = "MDD",
    "asd" = "ASD",
    "adhd_ipsych" = "ADHD",
    "insomn2" = "Insomnia",
    "eduAttainment" = "Education",
    "intel" = "IQ",
    "alcoholism_2019" = "Alcoholism",
    "ocd" = "OCD",
    "tourette" = "Tourettes",

    
    "alzBellenguezNoApoe" = "AD",
    "ms" = "MS",
    "pd_without_23andMe" = "PD",
    "epilepsyFocal" = "Epilepsy",
    "migraines_2021" = "Migraines",
    "stroke" = "Stroke",
    "als2021" = "ALS",
    
    
    "obesity" = "Obesity",
    "dm2" = "T2D",
    "lipidCholTotal" = "Cholesterol",
    "ra" = "RArthritis",
    "ibd" = "IBD",
    "uc" = "UC")

trait_category <- c(
    "sz3" = "Psychiatric",
    "bip2" = "Psychiatric",
    "mdd_ipsych" = "Psychiatric",
    "asd" = "Psychiatric",
    "adhd_ipsych" = "Psychiatric",
    "insomn2" = "Psychiatric",
    "eduAttainment" = "Psychiatric",
    "intel" = "Psychiatric",
    "alcoholism_2019" = "Psychiatric",
    "ocd" = "Psychiatric",
    "tourette" = "Psychiatric",

    
    "alzBellenguezNoApoe" = "Neurological",
    "ms" = "Neurological",
    "pd_without_23andMe" = "Neurological",
    "epilepsyFocal" = "Neurological",
    "migraines_2021" = "Neurological",
    "stroke" = "Neurological",
    "als2021" = "Neurological",
    
    
    "obesity" = "Other",
    "dm2" = "Other",
    "lipidCholTotal" = "Other",
    "ra" = "Other",
    "ibd" = "Other",
    "uc" = "Other")