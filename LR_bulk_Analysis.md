# Long-Read-RNAseq-Processing-Pipeline

This pipeline provides a fully containerized Singularity environment that bundles all required tools and dependencies for **long-read transcriptome analysis**. With a single command, the entire workflow from raw FASTQ input through sequence alignment, full-length isoform identification and quantification, and a rich set of downstream analyses can be executed reproducibly on any compatible system.

- **ONT (Oxford Nanopore)** — full-length cDNA and direct RNA reads, transcript isoforms are discovered and quantified with **bambu**.
- **PacBio (Iso-seq / HiFi)** — full-length reads are processed with the official **Iso-seq** suite.

Both branches converge on **minimap2** for sequence alignment, and share the same downstream modules (isoform classification, differential expression, differential transcript usage, alternative splicing, coding-potential prediction, fusion detection). **The pipeline now supports multiple samples and optional treatment-vs-control comparisons.**

---

# Part I Workflow

A schematic of the full long-read RNA-seq workflow. The pipeline branches at preprocessing according to `platform`, then merges at the alignment stage.

```
                                  ┌──────────────────────────┐
        raw reads (FASTQ / BAM)   │   Step 0  Raw read QC     │
                                  │   NanoPlot / NanoComp     │
                                  └────────────┬─────────────┘
                                               │
                 ┌─────────────────────────────┴─────────────────────────────┐
                 │                                                             │
        platform = "ONT"                                            platform = "PacBio"
                 │                                                             │
   ┌─────────────▼──────────────┐                            ┌────────────────▼─────────────────┐
   │ Step 1  Full-length reads  │                            │ Step 1  Full-length reads          │
   │   Pychopper                │                            │   lima  →  isoseq refine (→ FLNC)  │
   │   (orient + trim primers)  │                            │   isoseq cluster (→ HQ transcripts)│
   └─────────────┬──────────────┘                            └────────────────┬─────────────────┘
                 │                                                             │
                 └──────────────────────────┬──────────────────────────────────┘
                                             │
                            ┌────────────────▼─────────────────┐
                            │ Step 2  Splice-aware alignment    │
                            │   minimap2 -ax splice             │
                            │   samtools sort / index / flagstat│
                            └────────────────┬─────────────────┘
                                             │
              ┌──────────────────────────────┴──────────────────────────────┐
              │                                                              │
   ONT: Step 3  Transcript discovery & quant            PacBio: Step 3  Isoform collapse & quant
        bambu (novel isoforms + counts)                       isoseq collapse  →  pigeon (classify/filter)
              │                                                              │
              └──────────────────────────────┬──────────────────────────────┘
                                             │
                            ┌────────────────▼─────────────────┐
                            │ Step 4  Isoform classification    │
                            │   SQANTI3 (FSM/ISM/NIC/NNC, QC)   │
                            └────────────────┬─────────────────┘
                                             │
        ┌──────────────────────┬─────────────┼──────────────┬───────────────────────┐
        │                      │             │              │                       │
  Step 5 DE/DTE         Step 6 DTU /    Step 7 AS      Step 8 Coding         Step 9 Fusion / polyA
  DESeq2 / edgeR        switching       SUPPA2         potential             JAFFAL / tailfindr
                        IsoformSwitch                  CPAT / TransDecoder
        │                      │             │              │                       │
        └──────────────────────┴─────────────┼──────────────┴───────────────────────┘
                                             │
                            ┌────────────────▼─────────────────┐
                            │ Step 10  Sample QC & integration  │
                            │   PCA + correlation + HTML report │
                            └──────────────────────────────────┘
```

> The PNG/SVG version of this diagram can be added under `figures/` as in the End-seq pipeline.

---

# Part II Requirements

1.  **Recommended Specs**:
      * 16-core CPU (long-read alignment and Iso-seq clustering are CPU-heavy)
      * 64 GB RAM (bambu / Iso-seq cluster peak usage scales with read number)
      * Fast local storage; HiFi/ONT BAMs are large

2.  **Singularity**: Must be installed on your system. The installation steps are identical to the End-seq pipeline. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
            libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        cd /mnt/share/software
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1
        ./mconfig
        cd builddir
        make
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        singularity --version
        singularity -h
        ```

3.  **snakemake**: Snakemake must be installed and requires a Python 3 distribution.

      ```bash
      pip install snakemake
      ```

4.  **Reference Data**: Long-read transcript analysis needs (a) a genome FASTA, (b) a reference annotation GTF, and (c) primer FASTA files for the relevant platform. Below are steps for the human hg38 / GENCODE v46 reference. For other genomes, download the corresponding files and replace them as needed.

      ```bash
      mkdir basement_data
      cd basement_data

      # Genome FASTA + annotation GTF (GENCODE)
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz
      gunzip GRCh38.primary_assembly.genome.fa.gz
      gunzip gencode.v46.primary_assembly.annotation.gtf.gz

      # (Optional) keep only main chromosomes for genome FASTA
      awk '/^>/ {p=0} /^>chr[0-9XYM]/ {p=1} p' \
          GRCh38.primary_assembly.genome.fa > GRCh38.primary_assembly.genome.chr.fa
      samtools faidx GRCh38.primary_assembly.genome.chr.fa   # bambu/SQANTI need .fai

      # (Optional) pre-build minimap2 splice index to save time on repeated runs
      singularity exec --cleanenv EndSeq_or_LongRead.sif \
          minimap2 -x splice -d hg38.splice.mmi GRCh38.primary_assembly.genome.chr.fa
      ```

      For **CPAT** coding-potential prediction, also download the prebuilt human logit model and hexamer table from the CPAT data repository, and place them under `References/cpat/`.

5.  **Required Files**:

      ```bash
      project_directory/
      ├── Scripts/
      │     ├── scripts/                # generate_report.py, generate_figures.py,
      │     │                           # compile_html_report.py, run_bambu.R, run_isoswitch.R ...
      │     ├── config.yaml
      │     └── LongRead.smk
      ├── Containers/
      │     └── LongRead.sif
      ├── References/
      │     ├── GRCh38.primary_assembly.genome.chr.fa
      │     ├── GRCh38.primary_assembly.genome.chr.fa.fai
      │     ├── gencode.v46.primary_assembly.annotation.gtf
      │     ├── hg38.splice.mmi                 # optional pre-built minimap2 index
      │     ├── isoseq_primers.fasta            # PacBio Iso-seq 5'/3' primers (for lima)
      │     └── cpat/
      │           ├── Human_Hexamer.tsv
      │           └── Human_logitModel.RData
      ```

      - **LongRead.smk** — The main Snakemake workflow script.
      - **config.yaml** — Configuration file with paths, parameters, and sample information. ⚠️ Must sit in the same directory as `LongRead.smk`.
      - **scripts** — Helper scripts (report generation, figure generation, and R drivers for bambu / IsoformSwitchAnalyzeR / SUPPA2 post-processing). ⚠️ Must sit in the same directory as `LongRead.smk`.
      - **LongRead.sif** — Singularity container with all required software pre-installed: `NanoPlot`, `NanoComp`, `pychopper`, `minimap2`, `samtools`, `bambu (R)`, `isoseq`, `lima`, `pbmm2`, `pigeon`, `SQANTI3`, `gffread`, `DESeq2/edgeR (R)`, `DRIMSeq/IsoformSwitchAnalyzeR (R)`, `SUPPA2`, `CPAT`, `TransDecoder`, `JAFFAL`, `tailfindr (R)`, `deeptools`, `multiqc`.
      - **isoseq_primers.fasta** — PacBio Iso-seq primer sequences for `lima`; replace with your kit's primers if needed.
      - **cpat/** — CPAT model files for coding-potential prediction; swap species-specific models as needed.

---

# Part III Running

   * **Step 1: Edit `config.yaml`**

      The pipeline supports **multiple samples** and optional **condition comparisons**. Set `platform` to switch between ONT and PacBio preprocessing/quantification. All samples are listed in `sample_info`; comparisons are defined separately in `comparison_info` and may be omitted.

      ```yaml
      # ==========================================
      # Long-Read RNA-seq Pipeline Configuration
      # ==========================================

      # Sequencing platform: "ONT" or "PacBio".
      #   ONT    -> Pychopper + bambu
      #   PacBio -> Iso-seq (lima/refine/cluster/collapse/pigeon)
      platform: "ONT"

      # For ONT: "cDNA" or "dRNA" (direct RNA). Controls minimap2 / Pychopper options.
      # For PacBio this field is ignored.
      ont_library: "cDNA"

      # Pychopper cDNA kit (ONT only). e.g. PCS111, PCS110, LSK114.
      pychopper_kit: "PCS111"

      # List all samples. For ONT use FASTQ; for PacBio use the HiFi (CCS) BAM.
      sample_info:
        - sample: "Ctrl1"
          input: "/project_directory/rawdata/Ctrl1.fastq.gz"
          group: "Ctrl"
        - sample: "Ctrl2"
          input: "/project_directory/rawdata/Ctrl2.fastq.gz"
          group: "Ctrl"
        - sample: "Treat1"
          input: "/project_directory/rawdata/Treat1.fastq.gz"
          group: "Treat"
        - sample: "Treat2"
          input: "/project_directory/rawdata/Treat2.fastq.gz"
          group: "Treat"

      # Optional: condition-vs-condition pairs for DE / DTU / differential splicing.
      # Keys are comparison names; treatment/control reference 'group' values above.
      # Omit this section to skip differential analyses.
      comparison_info:
        Treat_vs_Ctrl:
          treatment: "Treat"
          control: "Ctrl"

      outputdir:    "/project_directory/result"
      genomeFa:     "/project_directory/References/GRCh38.primary_assembly.genome.chr.fa"
      annotationGtf:"/project_directory/References/gencode.v46.primary_assembly.annotation.gtf"
      minimapIndex: "/project_directory/References/hg38.splice.mmi"   # optional; else genomeFa is used
      isoseqPrimers:"/project_directory/References/isoseq_primers.fasta"   # PacBio only
      cpatDir:      "/project_directory/References/cpat"
      sif:          "/project_directory/Containers/LongRead.sif"

      threads: 16

      # Genome name for SQANTI/MACS-style parameter inference and report display
      # Supported: hg19/hg38 -> human, mm10/mm39 -> mouse, etc.
      genome_name: "hg38"

      # Which downstream modules to run (toggle on/off)
      modules:
        sqanti:         true     # isoform classification & QC
        de:             true     # gene/transcript differential expression (DESeq2/edgeR)
        dtu:            true     # differential transcript usage (IsoformSwitchAnalyzeR)
        splicing:       true     # alternative splicing events (SUPPA2)
        coding:         true     # coding-potential / ORF (CPAT + TransDecoder)
        fusion:         false    # fusion-transcript detection (JAFFAL)
        polya:          false    # poly-A tail length (ONT only: tailfindr)

      # Statistical thresholds
      de_padj:   0.05
      de_lfc:    1.0
      dtu_padj:  0.05

      # Report language: "en" or "cn"
      report_language: "en"
      ```

   * **Step 2: Run snakemake**

      ```bash
      snakemake -s LongRead.smk --use-singularity --cores 16 \
        --singularity-args "--bind /project_directory:/project_directory"
      ```

   * **Command Parameters**

      **edit `config.yaml`**
      - `platform`: `"ONT"` or `"PacBio"`. Selects the preprocessing/quantification branch. **(required)**
      - `ont_library`: `"cDNA"` or `"dRNA"`. ONT only; tunes minimap2/Pychopper. **(optional, default: cDNA)**
      - `pychopper_kit`: cDNA kit name for Pychopper. ONT only. **(optional)**
      - `sample_info`: List of all samples. Each item needs `sample`, `input` (FASTQ for ONT / HiFi BAM for PacBio), and `group`. **(required)**
      - `comparison_info`: Optional condition pairs for differential analyses. Omit to skip. **(optional)**
      - `outputdir`: Output directory. **(required)**
      - `genomeFa`: Genome FASTA (indexed `.fai`). **(required)**
      - `annotationGtf`: Reference annotation GTF. **(required)**
      - `minimapIndex`: Pre-built minimap2 `.mmi`; if absent, `genomeFa` is used directly. **(optional)**
      - `isoseqPrimers`: Iso-seq primer FASTA. PacBio only. **(required for PacBio)**
      - `cpatDir`: Directory with CPAT model files. **(required if `modules.coding: true`)**
      - `sif`: Path to the Singularity image. **(required)**
      - `threads`: Threads per rule. **(required)**
      - `genome_name`: Assembly name, used by SQANTI3 and shown in the report. **(required)**
      - `modules`: Per-module on/off switches for the downstream analyses. **(optional)**
      - `de_padj` / `de_lfc` / `dtu_padj`: Significance thresholds. **(optional)**
      - `report_language`: `"en"` or `"cn"`. **(optional, default: en)**

      **run snakemake** (identical meaning to the End-seq pipeline)
      - `--use-singularity`: Run each rule inside the container.
      - `--singularity-args`: Extra Singularity arguments (e.g. `--bind`, `--nv`).
      - `--cores`: Max parallel CPU cores.
      - `--bind`: Directories mounted into the container; include rawdata, scripts, container, and references. **(required)**

---

# Part IV Output

   * **Output Structure**

      ```bash
      result/
      ├── qc/
      │     ├── Ctrl1.NanoPlot/              # per-sample raw read QC (read length / quality)
      │     ├── NanoComp_report.html         # cross-sample comparison
      │     └── ...
      ├── preprocess/
      │     ├── Ctrl1.full_length.fastq.gz    # ONT: Pychopper full-length reads
      │     ├── Ctrl1.pychopper.stats         # ONT: classification stats
      │     ├── Ctrl1.flnc.bam                # PacBio: full-length non-concatemer reads
      │     ├── Ctrl1.refine.report.csv       # PacBio: isoseq refine report
      │     └── ...
      ├── bam/
      │     ├── Ctrl1.sorted.bam
      │     ├── Ctrl1.sorted.bam.bai
      │     ├── Ctrl1.flagstat.txt
      │     └── ...
      ├── transcripts/
      │     ├── ONT_bambu/                    # ONT branch
      │     │     ├── extended_annotations.gtf
      │     │     ├── counts_transcript.txt
      │     │     ├── counts_gene.txt
      │     │     └── CPM_transcript.txt
      │     └── PacBio_isoseq/                # PacBio branch
      │           ├── Ctrl1.clustered.bam
      │           ├── collapsed.gff
      │           ├── collapsed.abundance.txt
      │           └── ...
      ├── sqanti/
      │     ├── *_classification.txt          # FSM/ISM/NIC/NNC categories
      │     ├── *_corrected.gtf
      │     ├── SQANTI3_report.html
      │     └── *_filtered.gtf
      ├── differential/
      │     ├── DE/
      │     │     ├── Treat_vs_Ctrl.gene.DESeq2.tsv
      │     │     └── Treat_vs_Ctrl.transcript.DESeq2.tsv
      │     ├── DTU/
      │     │     ├── Treat_vs_Ctrl.isoformSwitch.tsv
      │     │     └── Treat_vs_Ctrl.switch_consequences.tsv
      │     └── splicing/
      │           ├── events.ioe              # SUPPA2 AS events
      │           └── Treat_vs_Ctrl.dpsi      # differential PSI
      ├── coding/
      │     ├── cpat.ORF_prob.tsv             # coding-potential scores
      │     └── transdecoder/                 # predicted ORFs / peptides
      ├── fusion/                             # if modules.fusion: true
      │     └── jaffa_results.csv
      ├── polya/                              # if modules.polya: true (ONT)
      │     └── Ctrl1.polya_lengths.tsv
      ├── figure/
      │     ├── transcript_PCA.pdf
      │     ├── transcript_cor.pdf
      │     ├── sqanti_category_barplot.pdf
      │     ├── isoform_per_gene.pdf
      │     └── ...
      ├── multiqc/
      │     ├── multiqc_data/
      │     └── multiqc_report.html
      ├── report/
      │     ├── LongRead_Analysis_Report.html
      │     ├── data/
      │     │     └── sample_statistics.json
      │     └── figures/
      └── log/
      ```

   * **Output Interpretation**

      Below each output is described with the **tool that produced it** and **how to use it**.

      ### Step 0 — Raw read QC (`NanoPlot` / `NanoComp`)
      - **`qc/<sample>.NanoPlot/`** — Per-sample read-length and read-quality distributions, N50, total yield. Long-read QC differs from short-read FastQC: the key metrics are read-length distribution and per-read quality (Qscore), not per-base quality.
      - **`qc/NanoComp_report.html`** — Side-by-side comparison of all samples' yield/length/quality. Use it to spot low-yield or degraded libraries before committing to alignment.

      ### Step 1 — Full-length read identification
      **ONT branch — `Pychopper`**
      - **`preprocess/<sample>.full_length.fastq.gz`** — Reads in which both 5' and 3' cDNA primers were found; reads are re-oriented to genomic sense and primers trimmed. Only full-length, correctly oriented reads pass to alignment. Command:
        ```bash
        pychopper -k PCS111 -t 16 -r report.pdf -S stats.tsv \
            input.fastq full_length.fastq
        ```
      - **`preprocess/<sample>.pychopper.stats`** — Counts of full-length / rescued / unclassified reads. A high full-length fraction indicates a good cDNA library.

      **PacBio branch — `lima` + `isoseq refine`**
      - **`preprocess/<sample>.flnc.bam`** — Full-Length Non-Concatemer reads: primers removed (`lima`), polyA trimmed and concatemers removed (`isoseq refine`). This is the canonical Iso-seq input for clustering. Commands:
        ```bash
        lima --isoseq --peek-guess -j 16 movie.hifi.bam isoseq_primers.fasta demux.bam
        isoseq refine --require-polya demux.5p--3p.bam isoseq_primers.fasta flnc.bam
        ```
      - **`preprocess/<sample>.refine.report.csv`** — Per-read polyA length and concatemer status.

      ### Step 2 — Splice-aware alignment (`minimap2` + `samtools`)
      - **`bam/<sample>.sorted.bam` / `.bai`** — Spliced genome alignment. minimap2 is run in splice mode, with platform-specific presets:
        ```bash
        # ONT cDNA (after Pychopper)
        minimap2 -ax splice -uf -k14 -t 16 ref.fa full_length.fastq | samtools sort -o sorted.bam
        # ONT direct RNA
        minimap2 -ax splice -uf -k14 ...
        # PacBio HiFi / Iso-seq
        minimap2 -ax splice:hq -uf -t 16 ref.fa hifi.fastq | samtools sort -o sorted.bam
        ```
      - **`bam/<sample>.flagstat.txt`** — `samtools flagstat` summary (total, mapped, supplementary, secondary). For long reads, watch the supplementary-alignment fraction (chimeric/split mappings) and the overall mapping rate.

      ### Step 3 — Transcript discovery & quantification
      **ONT branch — `bambu`**
      - **`transcripts/ONT_bambu/extended_annotations.gtf`** — Reference annotation augmented with novel isoforms discovered across all samples jointly.
      - **`transcripts/ONT_bambu/counts_transcript.txt` / `counts_gene.txt`** — Transcript- and gene-level count matrices (samples in columns). These feed the differential-expression and DTU modules. bambu is run from R:
        ```r
        library(bambu)
        ann <- prepareAnnotations("annotation.gtf")
        se  <- bambu(reads = bam_files, annotations = ann, genome = "genome.fa", ncore = 16)
        writeBambuOutput(se, path = "transcripts/ONT_bambu/")
        ```

      **PacBio branch — `isoseq cluster` / `collapse` + `pigeon`**
      - **`transcripts/PacBio_isoseq/collapsed.gff`** — Genome-mapped, redundancy-collapsed unique isoforms.
      - **`transcripts/PacBio_isoseq/collapsed.abundance.txt`** — Full-length read counts per collapsed isoform (the PacBio quantification). Commands:
        ```bash
        isoseq cluster2 flnc.bam clustered.bam
        pbmm2 align --preset ISOSEQ --sort ref.fa clustered.bam aligned.bam   # or minimap2 -ax splice:hq
        isoseq collapse aligned.bam collapsed.gff
        pigeon sort collapsed.gff -o collapsed.sorted.gff
        pigeon classify collapsed.sorted.gff annotation.sorted.gtf ref.fa
        pigeon filter collapsed.sorted_classification.txt
        ```
        `pigeon` is the SMRT-Link isoform classifier/filter (SQANTI-derived) — it removes likely artifacts (intra-priming, RT-switching).

      ### Step 4 — Isoform classification & QC (`SQANTI3`)
      - **`sqanti/*_classification.txt`** — Each isoform labelled relative to the reference: **FSM** (full splice match), **ISM** (incomplete), **NIC** (novel in catalog), **NNC** (novel not in catalog), antisense, intergenic, etc., plus QC flags (RT-switching, non-canonical junctions, intra-priming).
      - **`sqanti/SQANTI3_report.html`** — Visual QC report (junction support, length distributions, category fractions). This is the single most informative "is my transcriptome real?" report. Run:
        ```bash
        sqanti3_qc.py isoforms.gtf annotation.gtf genome.fa -d sqanti/ --report html
        sqanti3_filter.py rules sqanti/*_classification.txt --gtf isoforms.gtf
        ```
      - **`figure/sqanti_category_barplot.pdf`** — Aggregated category counts per sample (analogous to the End-seq peak-count barplot).

      ### Step 5 — Differential expression (`DESeq2` / `edgeR`)
      - **`differential/DE/*_vs_*.gene.DESeq2.tsv`** and **`*.transcript.DESeq2.tsv`** — log2 fold-change, p-value, and adjusted p-value per gene/transcript. Built from the bambu (ONT) or collapsed abundance (PacBio) count matrices. Only produced when `comparison_info` is set.

      ### Step 6 — Differential transcript usage / isoform switching (`IsoformSwitchAnalyzeR`, `DRIMSeq`)
      - **`differential/DTU/*_vs_*.isoformSwitch.tsv`** — Transcripts whose **relative usage** within a gene changes between conditions (even if total gene expression does not) — a phenomenon long reads are uniquely well-suited to detect.
      - **`differential/DTU/*_vs_*.switch_consequences.tsv`** — Predicted functional consequences of switches (ORF gain/loss, NMD sensitivity, domain loss), annotated by IsoformSwitchAnalyzeR.

      ### Step 7 — Alternative splicing events (`SUPPA2`)
      - **`differential/splicing/events.ioe`** — AS events derived from the (extended) annotation: SE (skipped exon), A5/A3 (alt 5'/3' splice site), MX (mutually exclusive), RI (retained intron), AF/AL (alt first/last exon).
      - **`differential/splicing/*_vs_*.dpsi`** — Differential PSI (ΔPSI) and significance per event. Commands:
        ```bash
        suppa.py generateEvents -i annotation.gtf -o events -e SE SS MX RI FL -f ioe
        suppa.py psiPerEvent --ioe-file events.ioe --expression-file tpm.tsv -o psi
        suppa.py diffSplice -m empirical -i events.ioe -p Ctrl.psi Treat.psi \
            -e Ctrl.tpm Treat.tpm -o Treat_vs_Ctrl
        ```

      ### Step 8 — Coding potential & ORF prediction (`CPAT` + `TransDecoder`)
      - **`coding/cpat.ORF_prob.tsv`** — Coding-probability score per transcript (CPAT), classifying novel isoforms as coding vs non-coding (potential lncRNA).
      - **`coding/transdecoder/`** — Predicted ORFs and peptide sequences (TransDecoder), useful for proteomic follow-up and for confirming productive vs NMD-target isoforms. The transcript FASTA is extracted from the GTF with `gffread`:
        ```bash
        gffread -w transcripts.fa -g genome.fa isoforms.gtf
        cpat.py -x Human_Hexamer.tsv -d Human_logitModel.RData -g transcripts.fa -o cpat
        TransDecoder.LongOrfs -t transcripts.fa && TransDecoder.Predict -t transcripts.fa
        ```

      ### Step 9 — Fusion & poly-A (optional)
      - **`fusion/jaffa_results.csv`** (`JAFFAL`) — Candidate fusion transcripts detected directly from long reads, with breakpoints and supporting read counts. Enable with `modules.fusion: true`.
      - **`polya/<sample>.polya_lengths.tsv`** (`tailfindr`/`nanopolish polya`, ONT only) — Per-read estimated poly-A tail length, a feature only accessible to long-read (especially direct-RNA) data. Enable with `modules.polya: true`.

      ### Step 10 — Sample-level QC & integrated report
      - **`figure/transcript_PCA.pdf`** — PCA of samples on the transcript-count matrix; detects outliers and batch effects (same role as `BW_compare_PCA.pdf` in End-seq).
      - **`figure/transcript_cor.pdf`** — Pairwise sample-correlation heatmap.
      - **`multiqc/multiqc_report.html`** — Aggregated NanoPlot/samtools/SQANTI metrics across all samples.
      - **`report/LongRead_Analysis_Report.html`** — Single-page integrated report (stat cards, summary table, compact figures), bilingual via `report_language`. Open in any browser for a full overview of the experiment.

---

# Reference

[1] Li, Heng. "Minimap2: pairwise alignment for nucleotide sequences." *Bioinformatics* 34.18 (2018): 3094–3100. doi:10.1093/bioinformatics/bty191

[2] Chen, Ying, et al. "Context-aware transcript quantification from long-read RNA-seq data with Bambu." *Nature Methods* 20 (2023): 1187–1195. doi:10.1038/s41592-023-01908-w

[3] PacBio. "Iso-Seq: Single-molecule full-length transcript sequencing." SMRT-Tools / pbbioconda. https://isoseq.how

[4] Pardo-Palacios, Francisco J., et al. "SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms." *Nature Methods* 21 (2024): 793–797. doi:10.1038/s41592-024-02229-2

[5] Vitting-Seerup, Kristoffer, and Albin Sandelin. "IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences." *Bioinformatics* 35.21 (2019): 4469–4471. doi:10.1093/bioinformatics/btz247

[6] Trincado, Juan L., et al. "SUPPA2: fast, accurate, and uncertainty-aware differential splicing analysis across multiple conditions." *Genome Biology* 19.1 (2018): 40. doi:10.1186/s13059-018-1417-1

[7] Davidson, Nadia M., et al. "JAFFAL: detecting fusion genes with long-read transcriptome sequencing." *Genome Biology* 23.1 (2022): 10. doi:10.1186/s13059-021-02588-5

