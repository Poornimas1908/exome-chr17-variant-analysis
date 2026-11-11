# Exome Variant Analysis: Chromosome 17 (GRCh38)



An end-to-end germline variant-calling workflow on a public exome from the Genome in a Bottle (GIAB) sample **HG003**, restricted to **chromosome 17 (GRCh38)**. The notebook takes raw reads through quality control, alignment, GATK best-practice processing, variant calling, filtering, population-frequency and functional annotation, and finishes by **benchmarking the calls against the GIAB truth set** and critically reviewing the candidate variants.



> \*\*Note:\*\* This is a learning and portfolio project on a public reference sample. It demonstrates a pipeline and its accuracy; it is not a clinical interpretation.

!\[Benchmark of raw and PASS calls against the GIAB truth set](results/figures/fig4\_benchmark.png)



## Contents

* [Highlights](#highlights)
* [Workflow](#workflow)
* [Results](#results)
* [Repository contents](#repository-contents)
* [Data sources](#data-sources)
* [How to run](#how-to-run)
* [Design decisions](#design-decisions)
* [Limitations](#limitations)
* [Possible extensions](#possible-extensions)
* [Acknowledgements](#acknowledgements)



## Highlights

|                                        **Metric**|                                           **Result**|
|-|-|
|Input|2,387,574 read pairs (151 bp, Illumina NovaSeq, IDT exome capture)|
|Reads kept after trimming|99.1% (Q30 94.0% → 94.4%)|
|Alignment|99.98% mapped, 99.5% properly paired|
|Duplicate rate|15.7%|
|Mean depth on capture targets|117x|
|Target bases at ≥20x|99.7% (48 of 11,933 regions below 20x)|
|Raw variant records|4,711 (4,060 SNVs, 647 indels, 4 mixed multiallelic sites)|
|Variants passing filters|4,381 (3,790 SNVs, 591 indels), Ti/Tv 2.66|
|Accuracy vs GIAB, PASS SNVs|precision 0.989, recall 0.976|
|Accuracy vs GIAB, PASS indels|precision 0.833, recall 1.000 (35 true-positive indels)|
|Rare or novel protein-changing candidates|51 (9 clustered, 4 with unusually high depth)|
|Reproducibility|Identical key counts in three independent runs (24 Oct, 31 Oct, 5 Nov 2026)|

## Workflow

BAM subset (chr17) → FASTQ → FastQC → fastp → BWA-MEM → SAMtools
   → GATK MarkDuplicates → BQSR → mosdepth (coverage)
   → GATK HaplotypeCaller → hard filters → bcftools norm
   → gnomAD allele frequencies → bcftools csq (functional consequences)
   → candidate selection → GIAB benchmark


|**Step**|**Tool**|**Purpose**|
|-|-|-|
|Read QC|FastQC, MultiQC|Check base quality, adapters and overall read health|
|Trimming|fastp|Remove adapters and low-quality reads (minimum length 36 bp)|
|Alignment|BWA-MEM, SAMtools|Map reads to chr17; sort and index the BAM|
|Duplicates and BQSR|GATK4|Flag PCR duplicates; recalibrate base qualities using known sites (Mills indels, 1000G SNPs)|
|Coverage|mosdepth|Depth over the capture targets|
|Variant calling|GATK HaplotypeCaller|Call SNVs and small indels in the targets (100 bp padding)|
|Filtering|GATK VariantFiltration|SNVs: QD < 2, FS > 60, SOR > 3, MQ < 40, MQRankSum < -12.5, ReadPosRankSum < -8. Indels: QD < 2, FS > 200, SOR > 10, ReadPosRankSum < -20|
|Normalisation and frequencies|bcftools norm, annotate|Split multiallelic sites; add gnomAD allele frequency|
|Functional annotation|bcftools csq|Predict consequences on Ensembl GRCh38 transcripts (release 113)|
|Benchmark|bcftools isec|Compare with GIAB truth in captured, high-confidence regions (match by position and alleles)|

**Tool versions (final run):** BWA 0.7.19, SAMtools and bcftools 1.25, GATK 4.6.2.0, fastp 1.4.0, mosdepth 0.3.3, plus FastQC, MultiQC, bedtools and Python (pandas, numpy, matplotlib). The first run (24 Oct) used SAMtools and bcftools 1.24 and mosdepth 0.3.14, with the same results. Exact versions are recorded in 'results/software\_versions.txt'.



## Results

### Read quality and alignment

fastp retained 99.1% of reads and raised the Q30 rate from 94.0% to 94.4%; 162,530 reads had adapters trimmed. After alignment, 99.98% of reads mapped and 99.5% were properly paired. The duplicate rate was 15.7%.

### Coverage

Mean depth over 2.16 Mb of capture targets (11,933 merged regions) was 117x. The fractions of target bases covered at ≥10x, ≥20x and ≥30x were 99.9%, 99.7% and 98.7%.

!\[Coverage over the capture targets](results/figures/fig2\_coverage.png)

### Variant calling and filtering

HaplotypeCaller produced 4,711 records: 4,060 SNVs, 647 indels and 4 mixed multiallelic sites (SNV plus indel, genotype 1/2) that 'SelectVariants' does not place in either the SNV or the indel set. Hard filtering retained 4,381 PASS variants (3,790 SNVs and 591 indels). The SNV transition/transversion ratio was 2.66, slightly below the roughly 2.8 to 3.0 usually reported for exomes, which is probably explained by the 100 bp padding (adding intronic and flanking sites) and a small number of remaining false positives. Heterozygous calls centre near an allele fraction of 0.5, as expected.

!\[Variant counts and allele balance](results/figures/fig3\_variants.png)

### Accuracy against the GIAB truth set (HG003, v4.2.1)

Evaluated only in regions that are both captured and GIAB-confident. Calls were matched to the truth set by position and alleles.

|**Callset**|**Type**|**TP**|**FP**|**FN**|**Precision**|**Recall**|
|-|-|-|-|-|-|-|
|Raw|SNV|1460|42|15|0.972|0.990|
|Raw|Indel|35|13|0|0.729|1.000|
|PASS|SNV|1439|16|36|0.989|0.976|
|PASS|Indel|35|7|0|0.833|1.000|

Filtering cut false-positive SNVs from 42 to 16 at a small cost in recall (0.990 → 0.976). The indel estimates rest on only 35 true-positive calls and are imprecise.

### Annotation

Of 4,434 normalised records, 3,994 (90%) are present in gnomAD. Taking the most severe consequence per variant (ignoring effects on NMD transcripts), 2,351 are intronic, 656 synonymous and 609 missense. High-impact classes are rare: 6 stop-gained, 9 splice donor or acceptor, 4 start-lost and 1 frameshift.

!\[Population frequency and consequences](results/figures/fig5\_annotation.png)

### Candidate variants

Of the 4,434 normalised records, 4,381 had depth ≥ 10 and genotype quality ≥ 20. Of these, 3,843 are common (gnomAD AF ≥ 1%), 123 are rare and 415 are absent from gnomAD. Keeping rare or novel variants with a protein-changing consequence gives **51 candidates** (41 missense, 6 in-frame insertions, 3 in-frame deletions, 1 frameshift; 50 distinct sites, because one site carries two insertion alleles).

Each candidate was flagged for two possible artifacts:

* **Clustered:** within 200 bp of another rare or novel variant (9 candidates)
* **High depth:** above 300x, roughly 2.5 times the mean target depth, which can indicate a collapsed repeat or duplicated region (4 candidates)

!\[Depth and quality of the candidate variants](results/figures/fig6\_candidates.png)

**Manual plausibility review.** Candidates were reviewed by position, genotype, depth, gnomAD frequency and local context. None has been validated (for example in IGV or ClinVar).

|Gene|Position (GRCh38)|Genotype|Depth|gnomAD AF|Assessment|
|-|-|-|-|-|-|
|FASN|chr17:82,088,256|0/1|204|0.000098|Credible: rare, not clustered, normal depth|
|USP43|chr17:9,728,628|0/1|240|0.000074|Credible: rare, not clustered, normal depth|
|AOC3|chr17:42,852,329|0/1|207|0.0067|Credible, but close to the 1% cutoff|
|SGSH|chr17:80,210,880|0/1|\~200|not in gnomAD|Credible: novel, not clustered; ClinVar check needed|
|CCDC40|chr17:80,090,216 and 80,090,265|1/1|343, 328|not in gnomAD|Likely mapping artifact: clustered, about 3x mean depth|
|KRTAP4-5|chr17:41,149,523|1/1|496|not in gnomAD|Likely artifact or unmatched polymorphism: repeat-rich gene, about 4x depth|
|POLR2A|chr17:7,513,774|1/1|152|not in gnomAD|Suspicious: a homozygous frameshift in an essential gene; inspect first|

The pattern is informative: the most dramatic candidates (homozygous, novel, very deep) are the least trustworthy, while the credible ones are less striking heterozygous calls. The full list is in `results/variants\_candidates.csv`.

## Repository contents

```
.
├── Exome\_chr17\_Variant\_Analysis.ipynb   # full analysis notebook (Google Colab)
├── results/
│   ├── variants\_all.csv                 # all 4,434 annotated records (consequence, gene, gnomAD AF)
│   ├── variants\_candidates.csv          # 51 rare/novel protein-changing candidates with flags
│   ├── benchmark\_vs\_giab.csv            # precision/recall table
│   ├── qc\_summary.csv                   # read and variant counts
│   ├── software\_versions.txt            # tool versions used
│   ├── figures/                         # figures used in the report (fig2 to fig6)
│   ├── compare\_raw/                     # bcftools isec output, raw calls vs truth
│   └── compare\_pass/                    # bcftools isec output, PASS calls vs truth
├── qc/
│   ├── multiqc\_report.html              # combined FastQC and fastp report
│   └── fastp.html                       # trimming report
└── README.md
```

Raw reads, BAM files and reference downloads are not included because of their size. The notebook downloads everything it needs.

## Data sources

|Resource|Source|
|-|-|
|Reads|HG003 NovaSeq exome (IDT capture, about 100x), from the DeepVariant exome case study test data; chr17 reads extracted and converted to FASTQ|
|Capture targets|IDT/Novogene GRCh38 capture BED (DeepVariant exome case study), chr17 regions merged|
|Reference|UCSC hg38, chromosome 17|
|BQSR known sites|Mills and 1000G gold-standard indels; 1000G phase 1 high-confidence SNPs (Broad hg38 bundle)|
|Allele frequencies|gnomAD ('af-only-gnomad.hg38', GATK best-practices resources)|
|Truth set|GIAB HG003 GRCh38 benchmark v4.2.1 and confident-region BED|
|Gene models|Ensembl GRCh38 release 113 (GFF3), chr17 converted to UCSC-style `chr17` naming|

## How to run

1. Open the notebook in Google Colab.
2. Run the first code cell. 'condacolab' restarts the runtime; then continue from the second code cell.
3. Run the remaining cells in order. A full run took about one hour on a standard Colab runtime.
4. Run the final export cell to zip the results, figures and QC reports.

Each pipeline stage has its own heading in the notebook, with a one-line explanation and the key result, so it can also be read without running it. Colab deletes its files when a session ends, so download the results straight after the run.

## Design decisions

* **Hard filters instead of VQSR.** A single sample is too small to train a recalibration model, so fixed thresholds recommended for small datasets were used.
* **'bcftools csq' instead of SnpEff.** SnpEff was planned, but its database could not be downloaded (the hosting server denied access in Colab and in a normal browser), so functional annotation uses 'bcftools csq' with Ensembl GFF3 gene models. It reports consequences on transcripts (coding, splice, intronic, UTR and non-coding) but not upstream, downstream or intergenic effects, and it gives no HIGH/MODERATE/LOW impact labels.
* **One consequence per variant.** 'csq' reports an effect for every overlapping transcript, so the most severe effect per variant (ignoring NMD transcripts) is used for counting and candidate selection.
* **Benchmark restricted to confident regions.** Truth-set comparison uses the intersection of the capture targets and the GIAB high-confidence regions, so errors in difficult regions do not distort the accuracy estimate.
* **Verification checks.** Variant counts were reconciled at each stage; for example, the 4 mixed multiallelic sites were identified with 'bcftools isec' to explain why 4,060 + 647 does not equal 4,711.

## Limitations

* **Single sample, one chromosome.** HG003 is a public reference individual, and the analysis covers chr17 only.
* **Hard filters** leave some false positives and false negatives.
* **Four mixed multiallelic sites** were not carried through filtering.
* **Indels in repeat regions are less reliable**, and the indel accuracy estimate rests on few calls.
* **"Not in gnomAD" can be overstated for indels**, because the same indel can be written in different ways and fail to match the database record.
* **Benchmarking method.** 'bcftools isec' matches records by position and alleles, does not compare genotypes and is not haplotype-aware. Dedicated tools such as hap.py or RTG vcfeval are recommended for benchmarking and may give somewhat different values, especially for indels and complex variants.
* **Benchmark scope.** Accuracy is measured only in captured, GIAB-confident regions.
* **Candidates are not validated.** Several apparent candidates show signs of mapping or repeat artifacts; none has been checked in IGV or ClinVar.
* **Functional annotation** is limited to 'bcftools csq' consequences (see Design decisions).

## Possible extensions

* Run the full exome (all chromosomes) or additional GIAB samples
* Compare hard filters with VQSR or a deep-learning caller such as DeepVariant
* Benchmark with hap.py or RTG vcfeval to obtain haplotype-aware, genotype-aware metrics
* Add SnpEff or VEP annotation for impact categories and ClinVar significance
* Inspect candidate variants in IGV and add screenshots
* Pass mixed multiallelic sites through filtering

## Acknowledgements

This project uses public resources from the Genome in a Bottle Consortium (GIAB), the DeepVariant exome case study, the Broad Institute GATK resource bundle, gnomAD, UCSC and Ensembl, and open-source tools developed by their respective authors.

## Author

Poornima S | www.linkedin.com/in/poornima-s-265596343

