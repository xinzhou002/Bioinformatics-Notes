# RStudio
# 批量解析 Sanger 测序 .ab1 文件
# 输出：每条序列的主链、次链序列，以及杂合位点列表

# 0. 环境准备 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("sangerseqR", quietly = TRUE))
  BiocManager::install("sangerseqR")

library(sangerseqR)

# 1. 单文件解析函数
process_ab1_file <- function(file_path) {
  tryCatch({
    # 读入 .ab1
    seq  <- readsangerseq(file_path)

    # 双链 base-call（ratio 可调节阈值）
    bc   <- makeBaseCalls(seq, ratio = 0.33)

    # 主链 / 次链序列
    primary   <- as.character(primarySeq(bc))
    secondary <- as.character(secondarySeq(bc))

    # 找杂合位点（主链字符 != 次链字符）
    het_idx <- which(strsplit(primary, "")[[1]] !=
                     strsplit(secondary, "")[[1]])
    het_txt <- if (length(het_idx)) paste(het_idx, collapse = ";") else "None"

    data.frame(
      FileName        = basename(file_path),
      PrimarySequence = primary,
      SecondarySequence = secondary,
      HetSites        = het_txt,
      stringsAsFactors = FALSE
    )

  }, error = function(e) {
    warning(sprintf("文件处理失败: %s\n%s", basename(file_path), e$message))
    data.frame(
      FileName        = basename(file_path),
      PrimarySequence = NA_character_,
      SecondarySequence = NA_character_,
      HetSites        = "ERROR",
      stringsAsFactors = FALSE
    )
  })
}

# 2. 批量处理函数
process_ab1_files <- function(directory_path) {
  ab1_files <- list.files(directory_path,
                          pattern = "\\.ab1$",
                          full.names = TRUE,
                          ignore.case = TRUE)

  if (length(ab1_files) == 0)
    stop("目录中未找到 .ab1 文件，请检查路径。")

  # lapply + rbind 组合为数据框
  res <- do.call(rbind, lapply(ab1_files, process_ab1_file))
  rownames(res) <- NULL
  res
}

# 3. 运行示例
directory_path <- "D:/job/test"   # ← 改成你的目录！

results <- process_ab1_files(directory_path)

# 4. 导出结果
output_path <- file.path(directory_path, "SangerSeqR_results.csv")
write.csv(results, output_path, row.names = FALSE)
cat("处理完成！结果已保存到:", normalizePath(output_path), "\n")
