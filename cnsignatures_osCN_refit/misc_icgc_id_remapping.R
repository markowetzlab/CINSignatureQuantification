library(dplyr)

# Load all data
icgc_segs <- readRDS("data/PCAWG_SNP6_segments_penalty70.rds")
icgc_samples <- read.table("data/sample.tsv",header = T,sep="\t")
pcawg_mapping <- read.table("data/PCAWG_summary_table_combined_annotations_v4.txt",header = T,sep = "\t")
activities <- readRDS("data/PCAWG_signature_activities_THRESH095_NAMESAPRIL21.rds")

# check activity matrix rownames against pcawg samplenames
activityIds <- rownames(activities)
all(activityIds %in% pcawg_mapping$samplename)
dim(pcawg_mapping[match(rownames(activities),pcawg_mapping$samplename),])

# check pcawg icgc_sample_ids against icgc_samples table
pcawg_mapping_icgc_ids <- pcawg_mapping_acts$icgc_sample_id
all(pcawg_mapping_icgc_ids %in% icgc_samples$icgc_sample_id)
dim(icgc_samples[match(pcawg_mapping_icgc_ids,icgc_samples$icgc_sample_id),])

# check segtable icgc_ids against icgc_sample_id
icgc_seg_ids <- unique(icgc_segs$sample)
all(icgc_seg_ids %in% icgc_samples$icgc_sample_id)
dim(icgc_samples[match(icgc_seg_ids,icgc_samples$icgc_sample_id),])

# left join tables pcawg and icgc_samples table on icgc_sample_id
merged_tables <- icgc_samples %>%
    full_join(.,pcawg_mapping,"icgc_sample_id")

dim(merged_tables)

# get samplename id from pcawg_mapping using mapping from icgc_sample_id
# check unique id length matches
length(merged_tables$samplename[match(icgc_seg_ids,merged_tables$icgc_sample_id)]) == length(unique(icgc_segs$sample))
# Confirm no NA mappings
any(is.na(merged_tables$samplename[match(icgc_seg_ids,merged_tables$icgc_sample_id)]))

# Add new sample col to segment table reflecting ids in activites matrix
icgc_segs$sampleNew <- merged_tables$samplename[match(icgc_segs$sample,merged_tables$icgc_sample_id)]

missingIds <- rownames(activities)[which(!rownames(activities) %in% unique(icgc_segs$sampleNew))]

all(!icgc_segs$sample %in% merged_tables$icgc_sample_id[merged_tables$samplename %in% missingIds])

length(unique(icgc_segs$sampleNew))
saveRDS(icgc_segs,"data/PCAWG_SNP6_segments_penalty70_remappedIDs.rds")
