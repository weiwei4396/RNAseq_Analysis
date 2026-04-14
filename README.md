事实上，目前所有的分析脚本直接询问Claude Code就可以基本解决了

Short-read and Long-read RNA-seq analysis
## Table of contents
- [LongReads](#LongReads)
- [ShortReads](#ShortReads)


## LongReads
### 1. Long reads Transcript Identification and Quantification
#### 1.1 [IsoQuant](https://github.com/ablab/IsoQuant) and [Spl-IsoQuant](https://github.com/algbio/spl-IsoQuant)
前者用来做bulk数据的识别和定量(这个算法也可以做单细胞或者空转, 但是不能区分umi, 如果已经做了barcode和umi的校正, 可以输入bam到这个算法中), 后者用来做单细胞和空转的识别和定量(只能输入fastq.gz, 我觉得他的barcode和umi识别和校正的效果不好)
<details>
<summary> </summary>

```python
isoquant.py \
  --threads 8 \
  --data_type nanopore \
  --reference /data/workdir/panw/C_code/sc_mouse_brain_data/Annotations_Fasta_Files/GRCm39.genome.fa \
  --genedb /data/workdir/panw/C_code/sc_mouse_brain_data/Annotations_Fasta_Files/mouse_gencode_vM34_comprehensive.gtf \
  --bam /data/workdir/panw/MouseHIPP3/P21/P21_M1_HIPP/P21_SC_M1_HIPP_BC_U8_sorted.bam \
  --read_group tag:BC \
  -o /data/workdir/panw/IsoQuant
```

</details>

#### 1.2 [Sicelore-2.1](https://github.com/ucagenomix/sicelore-2.1)
用来做长读长单细胞数据的定量
<details>
<summary> </summary>

- **1.扫描Nanopore的reads, 分配cell barcode;**

这一步扫描文库生成的cDNA序列, 对Nanopore fastq reads扫描polyA(T)和adaptor, 将其归类为pass; 扫描polyA(T)的长度参数大于等于15nt, A(T)≥75%, 在两端100nt以内, 找到polyA之后在下游找adaptor; 在得到第一次的pass的reads后, 使用cellranger生成的barcode list来分配barcode, 在给定barcode list的条件下, 99%的都会分配正确, 有正确barcode信息的归类于pass, 其余的为fail.

```shell
LongRead_fastqd=/data/workdir/panw/Data/LongBench/All_fastq/LR_sc_cellMix/SC_ONT.fastq.gz
bcumifinder=/data/workdir/panw/softwares/sicelore-2.1/Jar/NanoporeBC_UMI_finder-2.1.jar
outputdir=/data/workdir/panw/Data/LongBench/All_result/LR_ONT_GEX/Sicelore
cellRangerbarcodes=/data/workdir/panw/Data/LongBench/All_result/SC_GEX_results/outs/filter_barcodes_clean.tsv
java -jar -Xmx75g $bcumifinder scanfastq -d $LongRead_fastqd -o $outputdir --bcEditDistance 1 --cellRangerBCs $cellRangerbarcodes
```
- -jar参数表示Java的库应用程序的压缩文件;
- -Xmx参数表示为Java分配的最大内存大小; 50g表示分配50G, 500m表示500M;
- --bcEditDistance参数表示与barcode的距离;
- -d表示输入的fastq文件;
- -o表示输出文件保存的路径;
- --cellRangerBCs这个是可选的参数, 如果有barcode.tsv文件可以提供;
- 这段代码会自动占用服务器中所有的线程; 接近98个线程, 运行了一小时; 62G的fastq.gz文件;

- **2.Mapping**
这一步表示将pass_fastq文件使用minimap2比对到reference genome;
```shell
reference=/data/workdir/panw/Reference/hg38.fa
juncBED=/data/workdir/panw/Data/LongBench/All_result/LR_ONT_GEX/junctions.bed
tNum=30
minimap2 -ax splice -ub -k14 -w 4 --junc-bed $juncBED --sam-hit-only --secondary=no -t $tNum $reference $outputdir/passed/SC_ONT_passed.fastq.gz > $outputdir/SC_ONT.sam
samtools view -bS -@ $tNum $outputdir/SC_ONT.sam > $outputdir/SC_ONT.bam
samtools sort -@ $tNum -o ${outputdir}/SC_ONT_sorted.bam $outputdir/SC_ONT.bam
samtools index -@ $tNum ${outputdir}/SC_ONT_sorted.bam $outputdir/SC_ONT_sorted_index
```

- **3.UMI assignment**
在Nanopore bam文件中寻找UMI;
```shell
gtf=/data/workdir/panw/Reference/gencode_new.v40.gtf
java -jar -Xmx80g $bcumifinder assignumis --inFileNanopore ${outputdir}/SC_ONT_sorted.bam --outfile ${outputdir}/SC_ONT_sorted_umi.bam --ONTgene GE --annotationFile $gtf
```

- **4.Generate cell or spatial barcode x Gene-/Isoform-/Junction-level matrices**
```shell
sicelorejar=/data/workdir/panw/softwares/sicelore-2.1/Jar/Sicelore-2.1.jar
java -jar -Xmx75g $sicelorejar SelectValidCellBarcode I=BarcodesAssigned.tsv O=ValidBarcodes.csv MINUMI=1 ED0ED1RATIO=1
gtfToGenePred -genePredExt -geneNameAsName2 $gtf gencode.v38.refflat.txt
paste <(cut -f 12 gencode.v38.refflat.txt) <(cut -f 1-10 gencode.v38.refflat.txt) > gencode.v38.refFlat
java -jar -Xmx180g $sicelorejar IsoformMatrix I=${outputdir}/SC_ONT_sorted_umi.bam GENETAG=GE UMITAG=U8 CELLTAG=BC REFFLAT=gencode.v38.refFlat CSV=ValidBarcodes.csv DELTA=2 MAXCLIP=150 METHOD=STRICT AMBIGUOUS_ASSIGN=false OUTDIR=./UMI PREFIX=sicelore
```
</details>




## 2. SUPPA Alternative splicing events
只需要给定更新后的gtf(可以用BroCOLI识别后更新的GTF, 或者初始给定的GTF), [SUPPA](https://github.com/comprna/SUPPA) 就可以从gtf中得到可变剪切事件;

<details>
<summary> </summary>

```shell
python /data/workdir/panw/softwares/SUPPA2.3/suppa.py generateEvents -i updated_annotations.gtf -o HIPP_MX -f ioe -e MX --pool-genes
```
- i: 输入的gtf文件, 从这个文件中找可变剪切事件;
- o: 输出文件的名称;
- f: ioi表示转录本的可变剪切事件; ioe表示局部可变剪切事件, ioe才可以定义e参数;
- e: 表示具体类型的可变剪切事件, 如MX表示互斥外显子; SE, RI, FL(AF,AL), SS(A3,A5); 
- --pool-genes：这个参数可以充分利用novel的转录本, 因为novel转录本的gene_id命名并不是十分精确, 有些novel转录本应该与注释的是相同基因, 可能是由于开头和结尾部分超过基因范围导致命名写为NA, 而这个参数可以根据位置更精确地确定可变剪切事件;

</details>




## ShortReads
### 1. Transcriptome gene quantification and differential analysis
这个流程主要基于STAR做比对，然后用featurecount对基因表达做定量，最后使用DESeq2做差异分析。

(1) STAR对参考基因组创建索引
<details>
<summary> </summary>

```shell
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir star252b --genomeFastaFiles /data/workdir/reference/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa --sjdbGTFfile /data/workdir/reference/human/refdata-gex-GRCh38-2024-A/genes/genes.gtf --sjdbOverhang 149 --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 3
```
</details>

需要注意的是, 创建的基因组的索引在使用时需要与STAR的版本保持一致, 否则会报错! 上面代码中后面的4个参数与cellranger下载下来的reference使用的参数是一样的;

(2) 使用STAR比对
<details>
<summary> </summary>
  
```shell
for f1 in ${DataPath}/*_R1_001.fastq.gz
do
    sample=$(basename "${f1}" _R1_001.fastq.gz)
    f2="${DataPath}/${sample}_R2_001.fastq.gz"
    if [[ ! -f "${f2}" ]]; then
        echo "R2 file not found for ${sample}, skip."
        continue
    fi
    echo "======================================"
    echo "Processing sample: ${sample}"
    echo "R1: ${f1}"
    echo "R2: ${f2}"
    echo "======================================"

    $STAR \
      --runThreadN ${threads} \
      --runMode alignReads \
      --genomeDir "${STARgenome}" \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate \
      --sjdbOverhang 149 \
      --readFilesIn "${f1}" "${f2}" \
      --readFilesCommand zcat \
      --outFileNamePrefix "${STARresult}/${sample}." \
      --genomeSAindexNbases 14 \
      --outFilterMultimapNmax 10 \
      --outFilterMismatchNmax 10
done

```

</details>

--sjdbOverhang参数需要保证与创建索引时一致;

(3) featurecount对基因定量
<details>
<summary> </summary>
  
```shell
FeatureCount=/data/workdir/run_pengmx/subread-2.1.1-Linux-x86_64/bin/featureCounts
out_dir=/data/workdir/run_pengmx/FeatureCountResult
bam_files=$(find "$STARresult" -maxdepth 2 -type f -name "*.Aligned.sortedByCoord.out.bam" | sort)

echo "These bam files will be counted:"
echo "$bam_files"

"$FeatureCount" \
  -T 30 \
  -p \
  --countReadPairs \
  -t exon \
  -g gene_name \
  -a "$genomeGtf" \
  -o "$out_dir/pengmx_featureCounts_count.txt" \
  $bam_files

```

</details>

(4) R脚本DESeq2差异基因表达分析
<details>
<summary> </summary>
  
```R

# ============================================================
# DESeq2 差异表达分析流程
# 输入：featureCounts 输出文件 pengmx_featureCounts_count.txt
# ============================================================
# ---- 1. 加载包 ----
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


# ---- 2. 读入 featureCounts 结果 ----
# 第一行是 # Program 注释，所以 skip=1
raw <- read.table("D:/AAA博士汇报/pmx/pengmx_featureCounts_count.txt",
                  header = TRUE, sep = "\t",
                  skip = 1, check.names = FALSE)

# 前 6 列为注释：Geneid Chr Start End Strand Length
# 第 7 列开始是各样本 count
count_matrix <- raw[, 7:ncol(raw)]
rownames(count_matrix) <- raw$Geneid

# featureCounts 默认样本名是完整 BAM 路径，清理一下
colnames(count_matrix) <- gsub(".*/|\\.bam$", "", colnames(count_matrix))
head(count_matrix)
# 46-3:CRC46_3
colnames(count_matrix) <- c("CRC46_3", "CRC46_4", "CRC46_5", "CRC46_6", "CRC60_5", "CRC60_6",
                            "control_1", "control_2", "control_3", "control_4", "control_6")

# ---- 3. 构建样本信息表（coldata） ----
# 根据你实际的分组情况修改!
sample_names <- colnames(count_matrix)
sample_names   # 先打印出来看一眼

# 比如 CRC 样本名里都带 "CRC" 或 "T"（tumor）,正常样本带 "N" 或 "control"
# 根据实际命名规则来匹配,下面只是示例:
group <- ifelse(grepl("CRC|tumor|T$", sample_names), "CRC", "control")

coldata <- data.frame(
  row.names = sample_names,
  condition = factor(group, levels = c("control", "CRC"))
)
coldata   # 打印出来肉眼核对每个样本的分组对不对
# 一致性检查
stopifnot(all(rownames(coldata) == colnames(count_matrix)))
table(coldata$condition)   # 应该显示 control: 5, CRC: 6

# ---- 4. 构建 DESeq2 对象 ----
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData   = coldata,
                              design    = ~ condition)

# 过滤掉低表达基因（至少在若干样本里有 ≥10 reads）
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]

# ---- 5. 运行 DESeq2 ----
dds <- DESeq(dds)
resultsNames(dds)   # 查看对比项名字

# ---- 6. 提取结果 ----
res <- results(dds, contrast = c("condition", "CRC", "control"),
               alpha = 0.05)
summary(res)

# LFC 收缩(可视化更稳定)
# 先看一下 resultsNames(dds) 输出的名字,一般会是 "condition_CRC_vs_control"
library(DESeq2)
library(apeglm)
res_shrink <- lfcShrink(dds, coef = "condition_CRC_vs_control",
                        type = "apeglm")

# 排序并导出
res_df <- as.data.frame(res_shrink)
res_df$Geneid <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, "c:/Users/Lenovo/Desktop/DESeq2_CRC_vs_control_all.csv", row.names = FALSE)


# ---- 7. 筛选差异基因 ----
# 常用阈值：|log2FC| > 1 且 padj < 0.05
lgFd = 0.6
deg <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > lgFd)
deg$regulation <- ifelse(deg$log2FoldChange > 0, "Up", "Down")
table(deg$regulation)
write.csv(deg, "c:/Users/Lenovo/Desktop/DEGs_padj0.05_log2FC1.csv", row.names = FALSE)

# ---- 8. 可视化 ----
## 8.1 PCA（需要先做 vst/rlog 变换）
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition") +
  geom_text(aes(label = name), vjust = 2, size = 3) +
  theme_bw()
ggsave("c:/Users/Lenovo/Desktop/PCA.pdf", width = 6, height = 5)

## 8.2 样本间距离热图
sampleDists <- dist(t(assay(vsd)))
pheatmap(as.matrix(sampleDists),
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         filename = "c:/Users/Lenovo/Desktop/sample_distance_heatmap.pdf")

## 8.3 火山图
res_df$sig <- "NS"
res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange > lgFd] <- "Up"
res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange < -lgFd] <- "Down"

ggplot(res_df, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
  geom_vline(xintercept = c(-lgFd, lgFd), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(x = "log2 Fold Change", y = "-log10(padj)", title = "Volcano plot")
ggsave("c:/Users/Lenovo/Desktop/volcano.pdf", width = 6, height = 5)

## 8.4 Top 差异基因表达热图
top_genes <- head(rownames(deg)[order(deg$padj)], 30)
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)   # 中心化
pheatmap(mat,
         annotation_col = coldata,
         show_rownames = TRUE,
         cluster_cols = TRUE,
         filename = "c:/Users/Lenovo/Desktop/top30_DEG_heatmap.pdf")

```

</details>




















