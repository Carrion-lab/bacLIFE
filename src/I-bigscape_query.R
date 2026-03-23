library(DBI)
library(RSQLite)
library(data.table)

db_path <- "intermediate_files/BiG-SCAPE/bigscape_output/bigscape_output.db"

message("Querying MIBIG SQLiteDB")

# Connect to the database
con <- dbConnect(RSQLite::SQLite(), db_path)

# Define the query
query <- "
SELECT
    br.id,
    br.gbk_id,
    br.category,
    br.product,
    brf.family_id,
    f.bin_label,
    g.path,
    g.organism,
    g.description
FROM bgc_record br
LEFT JOIN bgc_record_family brf
    ON br.id = brf.record_id
LEFT JOIN Family f
    ON brf.family_id = f.id
LEFT JOIN gbk g
    ON br.gbk_id = g.id;
"

# Run the query and convert to data.table
result_dt <- as.data.table(dbGetQuery(con, query))

# Disconnect from the database
dbDisconnect(con)

message("Creating and writing GCF tables")

# Make a copy
dt2 <- copy(result_dt)

# Keep only rows where bin_label == "mix"
dt2 <- dt2[bin_label == "mix"]
dt2[, bin_label := NULL]

BGC_annot <- copy(dt2)

# Rename columns
setnames(dt2,
         old = c("category", "organism", "path", "family_id"),
         new = c("BGC Class", "Organism", "BGC Name", "GCF No"))

# Remove first two columns (id, gbk_id)
dt2 <- dt2[, -c(1, 2,4,8)]

# Reorder columns
setcolorder(dt2, c("BGC Name", "GCF No", "Organism", "BGC Class"))

#extract name
dt2[, `BGC Name` := sub(".*\\/([^\\/]+)\\.gbk$", "\\1", `BGC Name`)]
#Only annotations
annot <- unique(dt2[, .(`GCF No`, `BGC Class`)])

fwrite(dt2, "intermediate_files/BiG-SCAPE/GCF_annotation.txt", sep = "\t")
fwrite(annot, "intermediate_files/BiG-SCAPE/annotation.txt")

message("Preparing BGC annotations")

setnames(BGC_annot,
         old = c("category", "organism", "path", "family_id"),
         new = c("BGC Class", "Organism", "BGC Name", "GCF"), skip_absent=TRUE)

# Remove first two columns (id, gbk_id)
BGC_annot <- BGC_annot[, -c(1, 2)]

#extract name
BGC_annot[, `BGC Name` := sub(".*\\/([^\\/]+)\\.gbk$", "\\1", `BGC Name`)]

# Reorder columns
setcolorder(BGC_annot, c("GCF","BGC Name","Organism", "BGC Class", "product","description"))

#Rename GCF to match previous baclife version V.1.1 and below
BGC_annot$GCF <- paste0("GCF",BGC_annot$GCF)

write.table(BGC_annot, 'intermediate_files/BiG-SCAPE/BGC_descriptions.txt', col.names = T, row.names = F, quote = F, sep = '\t')

message("Bigscape DB processing done")
