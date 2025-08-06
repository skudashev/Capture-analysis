### Long-read Isoform Annotation Pipeline

This Snakemake pipeline processes long-read CaptureSeq data from raw basecalled FASTQs through isoform detection, read collapsing, and quality control. 

The pipeline uses a `samples.tsv` to control sample-specific input paths.

### 🛠 Tool Dependencies and Versions

This pipeline relies on pre-installed tools on the Earlham Institute cluster. Tools are sourced at the start of each batch job using `source package` or module-style commands.

| Tool | Version |
| --- | --- |
| **Snakemake** | 8.14.0 |
| **Minimap2** | 2.24 |
| **Samtools** | 1.15.1 |
| **SQANTI3** | 5.1 |
| **PicardTools** | 2.23.9 |
| **Restrander** | 1.0.0 |
| **Bedtools** | 2.30.0 |
| **SeqKit** | 2.5.1 |

### Environment Setup

I source my `~/.bashrc` which has conda, so you will need an env with pandas installed
