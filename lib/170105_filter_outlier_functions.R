filter_outlier <- function(input, output, barcode, deviation = 5){
  #[input] is the input count data, 2 columns are "seq" and "count"
  #[output] is the output count data, 2 columns are "seq" and "count"
  #[barcode] is the a data.frame that has three columns "seq", "rs", "nt"
  data <- merge(input, output, by.x = "seq", by.y = "seq")
  colnames(data) <- c("seq", "input_count", "output_count")
  
  data <- merge(data, barcode, by.x = "seq", by.y = "seq")
  data.split <- split(data, data[,c("rs","nt")], drop = TRUE)
  filter_outlier_for_each_rs <- function(count_table){
    #[count_table] is the count table for each rs, should have 1-2 alleles, 2 mostly, but also need to consider situations 
    change.mad <- mad(count_table$output_count/count_table$input_count)
    change.median <- median(count_table$output_count/count_table$input_count)
    lower <- change.median-deviation*change.mad
    upper <- change.median+deviation*change.mad
    count_table$lower <- lower
    count_table$upper <-upper
    count_table[(count_table$output_count/count_table$input_count) <= upper & (count_table$output_count/count_table$input_count) >= lower,]
  }
  data.split.filtered <- lapply(data.split, filter_outlier_for_each_rs)
  data.filtered <- do.call(rbind, data.split.filtered)
  data.filtered$input_CPM <- data.filtered$input_count/sum(data.filtered$input_count)*1000000
  data.filtered$output_CPM <- data.filtered$output_count/sum(data.filtered$output_count)*1000000
  data.filtered
}


