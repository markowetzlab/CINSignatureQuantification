icgc_segs <- readRDS("../../Downloads/PCAWG_SNP6_segments_penalty70.rds")
icgc_samples <- read.table("../../Downloads/sample.tsv",header = T,sep="\t")
pcawg_mapping <- read.table("../../Downloads/PCAWG_summary_table_combined_annotations_v4.txt",header = T,sep = "\t")
activities <- readRDS("../../Downloads/PCAWG_signature_activities_THRESH095_NAMESAPRIL21.rds")

activityIds <- rownames(activities)
pcawg_mapping_acts <- pcawg_mapping[pcawg_mapping$samplename %in% activityIds,]

pcawg_mapping_icgc_ids <- pcawg_mapping_acts$icgc_sample_id

icgc_samples_acts <- icgc_samples[icgc_samples$icgc_sample_id %in% pcawg_mapping_icgc_ids,]

icgc_samples_acts_filt <- icgc_samples_acts[icgc_samples_acts$project_code %in% c("OV-AU","ESAD-UK"),]

table(icgc_samples_acts_filt$project_code)
all(icgc_samples_acts_filt$icgc_sample_id %in% unique(icgc_segs$sample))
