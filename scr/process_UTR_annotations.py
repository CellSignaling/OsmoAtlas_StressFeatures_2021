import pandas as pd

# Read 3' UTR annotation skipping commented lines and the first line
utr3 = pd.read_table("data/derived_data/Nagalakshmi_2008_3UTRs_V64_wo_conflict_UTRs.bed", header = None, comment = '#', skiprows = 1)

# Eliminate chr from column 1 (chrI --> I)
utr3[0] = utr3[0].str.split("chr", expand=True)[1]
# Change mt label to Mito for compatibility with the genome
utr3[0] = utr3[0].replace("mt", "Mito")
# Change name column from YAL062W_3UTR to 3UTR:YAL062W
utr3[3] = utr3[3].str.split("_", expand = True)[1] + ":" + utr3[3].str.split("_", expand = True)[0]
# Write to a bed file
utr3.to_csv('data/derived_data/Nagalakshmi_2008_3UTRs_V64_formatted.bed', sep = '\t', index = False, header = False)

# Read 5' UTR annotation skipping commented lines and the first line
utr5 = pd.read_table("data/derived_data/Nagalakshmi_2008_5UTRs_V64_wo_conflict_UTRs.bed", header = None, comment = '#', skiprows = 1)
# Eliminate chr from column 1 (chrI --> I)
utr5[0] = utr5[0].str.split("chr", expand=True)[1]
# Change mt label to Mito for compatibility with the genome
utr5[0] = utr5[0].replace("mt", "Mito")
# Change name column from YAL062W_3UTR to 3UTR:YAL062W
utr5[3] = utr5[3].str.split("_", expand = True)[1] + ":" + utr5[3].str.split("_", expand = True)[0]
# Write to a bed file
utr5.to_csv('data/derived_data/Nagalakshmi_2008_5UTRs_V64_formatted.bed', sep = '\t', index = False, header = False)