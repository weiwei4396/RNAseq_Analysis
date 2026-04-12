## Transcriptome_Analysis
Long-read RNA-seq analysis
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




## 3. SUPPA Alternative splicing events
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





























