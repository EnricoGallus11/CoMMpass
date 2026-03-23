# !/usr/bin/r
# MMRF CoMMpass Expression Matrix — Load, Filter, Transpose
# Last update: 17/03/2026

# Adding Libraries (I've previously installed)
library(tidyverse)
library(data.table)

# Define the path to the directory containing the DataSet
input_dir <- "/Users/biocompute/Alma Mater Studiorum Università di Bologna/Bioinformatics Seràgnoli - expression_estimates_transcript_based"

# Read the TSV file using fread() and storing it into a variable
# paste0() concatenates the directory path with the filename (no need to write the whole path again)
input_df <- fread(paste0(input_dir, "/MMRF_CoMMpass_IA22_salmon_transcriptUnstrandedIgFiltered_counts.tsv"))

# str() shows the structure of the specified file
str(input_df)

# head() prints a part of the file to allow a quick inspection
head(input_df)

# View() opens the data as a window near the R script page
View(input_df)

processed_df <- input_df %>%
  as_tibble() %>%                           # Converting a table (.tsv) into a tibble makes it more manageable for R 
  distinct()                                # Keeps all the rows holding names unique (removes duplicates)
# And just like this we're now dealing with one single gene per row



# Now we want to check whether there are columns reporting unwanted infos
# To achieve so we base on the name of columns, as Viola told us this DF has an encoded pattern defining column content
# so we want to define a specific code that would allow us to see if there are column showing unwanted content
# the columname pattern is a whole string built as follows:
# "DiseaseCode"_"PatientNumber"_"TherapyCycle"_"SampleTissue"_"CloneUnderInvestgtn"
#       1              2               3              4                 5
# 1, 3, 4, 5 are meant to be kept constant, filtering out whatever differs
# 2 can change but should be kept unique after filtering above


colnames(processed_df)                      #  names of the columns

# ======================
# GREP(), ... OUR SAVIOR
# ======================
# grep allows looking for specific patterns within text, in the column name row in this case, how? grep(DataFrame)!!
# okay but how to define that pattern?
# There is a REALLY DAMN HARD encoding language to """help""" user scanning for complex patterns (as in this case):

#______________________
# BASIC COMBOS

# grep("text", text) --> scans for "text" within the variable text 
# this will return the position of that exact string section (looks like:   [8] like in lists indexing)

# grep("text", text, value=TRUE)
# will return THE FULL NAME OF THE COLUMNS containing "text" (seems useful)

# grep("text", text, invert=TRUE)
# returns whatever differs than specified :)

#________________________
# REGULAR EXPRESSIONS

# grep("^MMRF", DF, value=TRUE) 
# ^ <- will return all the strings in DF STARTING BY "MMRF" 

# grep("CD138pos$", DF, value=TRUE)
# $ <- will return all the strings in DF ENDING WITH "CD138pos"

# grep("hola", DF, value=TRUE, ignore.case=TRUE)
# will turn of case sensitive scanning and return al columns containing either hola, HOLA, HolA, ...

# grep("fontanesi_._", DF, value=TRUE)
# will return all columns with "fontanesi_WHICHEVERSINGLECHAR_"

# grep("fontanesi_[0-9]", DF, value=TRUE)
# will return all columns with "fontanesi_[ANYSINGLENUMBER]" :).      [0-9]+ = one or more digits

# grep("fontanesi_.*_", DF, value=TRUE)
# will return all columns with "fontanesi_WHICHEVERSEQUENCEOFCHARS_"

#_________________________

# SO I would trying use inverse to display all columns that are not showing BM, first cycle's patients nor CD138pos and see

grep("MMRF_[0-9][0-9][0-9][0-9]_1_BM_CD138pos", colnames(processed_df), invert=TRUE)

# okay so we have 168 columns that must be filtered out, but how many do we have for each case?

# Let's check how many are from multiple therapy cycles
grep("MMRF_[0-9][0-9][0-9][0-9]_1_.*_.*", colnames(processed_df), invert=TRUE)
# 154... 

# Let's check for the Sample's Tissue
grep("MMRF_[0-9][0-9][0-9][0-9]_[0-9]_BM_.*", colnames(processed_df), invert=TRUE)
# 19(-1)...

# Cell type's turn
grep("MMRF_[0-9][0-9][0-9][0-9]_[0-9]_.*_CD138pos", colnames(processed_df), invert=TRUE)
# 1 (the first one), so none, YAYYY

# Let's see the 19 Samples
View(data.frame(columns = colnames(processed_df)[grep("MMRF_[0-9][0-9][0-9][0-9]_[0-9]_BM_.*", colnames(processed_df), invert = TRUE)]))
# view a dataframe with column values corresponding to the ones from processed_df after grepping....
# So all these 19(-1) cases have PB (Peripheral Blood I guess) instead of BM

# At this point we just have to make all those columns disappear and transpose the matrix
# Instead of manually deleting all those columns (as I would do with Pandas in python)
# we'll initialize a new DF holding all the WANTED columns, so:

untrasp_df <- processed_df[, grep("MMRF_[0-9][0-9][0-9][0-9]_1_BM_CD138pos", colnames(processed_df), value=TRUE)]
head(untrasp_df)

# booom, let's do the same tests we've done above
grep("MMRF_[0-9][0-9][0-9][0-9]_1_.*_.*", colnames(untrasp_df), invert=TRUE)      #integer(0)
grep("MMRF_[0-9][0-9][0-9][0-9]_[0-9]_BM_.*", colnames(untrasp_df), invert=TRUE)  #integer(0)
grep("MMRF_[0-9][0-9][0-9][0-9]_[0-9]_.*_CD138pos", colnames(untrasp_df), invert=TRUE)  #integer(0)


# okay so now genes/transcripts are ROWS and patients are COLUMNS
# we want to flip diagonally the matrix so that rows are column and column are rows

# column_to_rownames() takes a column values and turns it into the row names 
# we'll use it on the 1st column (so the one containing all the column names)
# colnames(.)[1] just grabs whatever that first column is called 
# while t() is the transpose function but it doesn't return a tibble so we need to convert again like before

# as_tibble(rownames = "sample_id") does two things:
# 1) converts the matrix back into a tibble
# 2) takes the row names (which after transposing are the old column names = "DiseaseCode"_"PatientNumber"_"TherapyCycle"_"SampleTissue"_"CloneUnderInvestgtn")
#    and puts them into a proper column called "sample_id" instead of leaving them as undeclared names


transposed_df <- untrasp_df %>%
  column_to_rownames(var = colnames(.)[1]) %>%  #var is like saying take THIS and turn it into rownames
  t() %>%
  as_tibble(rownames = "sample_id")

# Error in `.rowNamesDF<-`(x, value = value) :   💕💕💕💕💕 Nice!
# duplicate 'row.names' are not allowed

# column_to_rownames() is complaining because the 1st column is NOT the gene IDs column anymore :(
# it's a patient sample column like MMRF_1021_1_BM_CD138pos 💀
# I deleted the column that was telling me what each row was 💀
# We need to manually add that first column back... (thankfully I created a new DF at each step)
# we save its name, grep the sample columns as before, then glue them together

# ==============
# c() FUNCTION:
# ==============
# Combine function creates a vector by combining specified values

# x <- c(1, 2, 3, 4)  will create an x vector storing those values [1, 2, 3, 4]
# Despite python, IF YOU MIX TYPES (like int, float, string) it forces everything to the most ""flexible"" type (string)

#Since we want to take values from one column we'll need something like: processed_df[, c(column_name)]

# We grab the name of that first column, before losing it (again), from processed_df and we store it into a new variable 
gene_id_col <- colnames(processed_df)[1]

# We do the same with grep in the wanted sample columns (exactly like before)
wanted_samples <- grep("MMRF_[0-9][0-9][0-9][0-9]_1_BM_CD138pos", colnames(processed_df), value=TRUE)

# c() concatenates: first the gene ID column name, then all the wanted sample names
# so now: COLUMN1 = gene IDs, COLUMN 23456+ = "DiseaseCode"_"PatientNumber"_"TherapyCycle"_"SampleTissue"_"CloneUnderInvestgtn")
untrasp_df <- processed_df[, c(gene_id_col, wanted_samples)]

# let's check: the first column should be "Transcript"
colnames(untrasp_df)[1]

head(untrasp_df)

# same sanity checks as before 
grep("MMRF_[0-9][0-9][0-9][0-9]_1_.*_.*", colnames(untrasp_df), invert=TRUE)
grep("MMRF_[0-9][0-9][0-9][0-9]_[0-9]_BM_.*", colnames(untrasp_df), invert=TRUE)
grep("MMRF_[0-9][0-9][0-9][0-9]_[0-9]_.*_CD138pos", colnames(untrasp_df), invert=TRUE) #(now they return 1 instead of integer(0), that 1 is the gene ID col)

# NOW let's transpose again
transposed_df <- untrasp_df %>%
  column_to_rownames(var = colnames(.)[1]) %>%    # THIS time the 1st col IS the gene IDs, so it works
  t() %>%                                          # flip the matrix diagonally
  as_tibble(rownames = "sample_id")                # patient IDs land in the "sample_id" column

View(transposed_df)
# VITTORIA