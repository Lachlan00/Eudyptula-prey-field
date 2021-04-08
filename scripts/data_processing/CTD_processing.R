# CTD cleaning - sorting messy data directory
# Moves all CTD cast data into clean directory and
# ensures no duplicate data
setwd("~/Development/PhD/repos/Eudyptula")

# input directory
input_dir <- 'data/CTD/sorted/all_files/'
output_dir <- 'data/CTD/processed/casts/'

# get the file, id and extension lists
file_ls <- list.files(path=input_dir, recursive=FALSE)
path_ls <- unlist(lapply(file_ls, function(x) paste0(input_dir, x)))
id_ls <- unlist(lapply(file_ls, function(x) str_split(x, '[.]')[[1]][1]))
ext_ls <- unlist(lapply(file_ls, function(x) str_split(x, '[.]')[[1]][2]))
dest_ls <- unlist(lapply(file_ls, function(x) paste0(output_dir, x)))

# create a data frame
df <- data.frame(file=file_ls,
                 id=id_ls,
                 ext=ext_ls,
                 path=path_ls,
                 dest=dest_ls)
# correct data formats
df$path <- as.character(df$path)
df$dest <- as.character(df$dest)

# check for duplicate cast ids
id_table <- table(factor(id_ls))
df$id_count <- as.integer(id_table[df$id])

# don't need to deal with xlsx or mat files so we can drop
df <- df[!df$ext %in% c('mat', 'xlsx'),]

# count casts
cast_no <- length(levels(df$id))
print(paste(cast_no, 'casts to process..'))

# move all files that are csv and id_count = 1
moveset <- df[df$id_count == 1 & df$ext == 'csv',]
df <- df[!(df$id_count == 1 & df$ext == 'csv'),]
file.copy(moveset$path, output_dir)
row.names(df) <- NULL

# find files with matching ids and copy across only the csv
# reset factors levels
df <- droplevels(df)
# check each id contains a csv file
df$id_csv <- sapply(df$id, function(x) 'csv' %in% df[df$id == x,]$ext)
# if csv present dump ctd and ctdx files
df <- df[!(df$id_csv & (df$ext %in% c('ctd', 'ctdx'))),]
# drop levels again
df <- droplevels(df)
# final check
if (length(levels(df$id)) == nrow(df)){
  print('All looks good!')
} else {
  print('Something may have gone wrong..')
}
# copy all remaining csv files
moveset <- df[df$ext == 'csv',]
df <- df[!df$ext == 'csv',]
file.copy(moveset$path, output_dir)
row.names(df) <- NULL
# confirm done
if (nrow(df) == 0){
  print('Success!')
} else {
  print('Something may have gone wrong..')
}
