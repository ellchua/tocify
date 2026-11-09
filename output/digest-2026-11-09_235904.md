# Weekly ToC Digest (week of 2026-11-09)

This digest is automatically generated from this week's RSS items and categorized into Neuroblastoma, AI, Methods, and Other.

**Included:** 34 (all ranked papers)  
**Scored:** 34 total items

**Models:** `gpt-5-nano` first-pass, `gpt-5.2` re-rank top 40

---

## Neuroblastoma (22 shown / 22 total)

### [PRISM: probabilistic cell-type resolution with integrated spatial modelling for simultaneous ambient RNA decontamination and doublet correction in spatial transcriptomics](https://www.biorxiv.org/content/10.64898/2026.11.07.765042v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.81**
Published: 2026-11-09
Tags: spatial transcriptomics, ambient RNA, doublet correction, EM algorithm, cell typing, QC

Introduces an EM framework that jointly handles two common spatial transcriptomics artifacts—ambient RNA contamination and segmentation doublets—while assigning probabilistic cell-type identities. This is highly transferable to spatial profiling of pediatric tumors and microenvironment analyses where these confounders can drive false signals.

<details>
<summary>RSS summary</summary>

Spatial transcriptomics count matrices are confounded by ambient RNA contamination and cell-segmentation doublets, two artifacts that co-occur but are rarely modelled jointly. We present PRISM, an Expectation-Maximization framework that simultaneously (i) assigns each cell a probabilistic cell-type identity, (ii) detects doublets via log-likelihood ratio testing of two-component mixtures, and (iii) estimates a per-cell contamination fraction anchored to a spatially resolved ambient profile. In a…

</details>

---

### [DeSCENT: Deconvolutional Single-Cell RNA-seq Enhances Transcriptome-based Cancer Survival Analysis](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag810/9016676?rss=1)
*Bioinformatics (Oxford Academic)*  
Score: **0.77**
Published: 2026-11-06
Tags: bulk RNA-seq, scRNA-seq, deconvolution, survival analysis, translational, biomarkers

Bridges the common gap between scRNA-seq (heterogeneity) and bulk cohorts (survival outcomes) by generating reference-informed cellular profiles from bulk RNA-seq for survival prediction. This is directly aligned with translational biomarker development using bulk clinical cohorts informed by single-cell atlases.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle">Motivation</div>Single-cell RNA-seq captures tumor heterogeneity, but it is rarely available together with patient survival outcomes. Bulk RNA-seq cohorts provide survival data but lack direct measurements of cellular expression.<div class="boxTitle">Results</div>DeSCENT generates reference-informed cellular profiles from bulk RNA-seq and combines both representations for survival prediction. We evaluated De…

</details>

---

### [THiRST: Computer Vision Empowered Three-dimensional High-Resolution Spatial Transcriptomics Analytics](https://www.biorxiv.org/content/10.64898/2026.10.29.762926v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.75**
Published: 2026-11-09
Tags: spatial transcriptomics, 3D reconstruction, computer vision, serial sections, tissue architecture

Presents a framework for reconstructing 3D spatial transcriptomic representations from serial sections using morphology-guided registration and spatial matching. This is relevant for understanding continuous tumor architecture and invasion niches in spatially profiled pediatric tumors.

<details>
<summary>RSS summary</summary>

High-resolution spatial transcriptomics (ST) preserves gene expression within tissue architecture but remains limited to individual two-dimensional (2D) sections, restricting analysis of continuous three-dimensional (3D) structures and cross-section consistency. We present THiRST, a computational framework for reconstructing multi-resolution 3D spatial transcriptomic representations from serial sections. THiRST combines coarse-to-fine morphology-guided registration, spatial matching of ST units,…

</details>

---

### [Structured multi-view relational learning for drug response prediction](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag619/9022080?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.74**
Published: 2026-11-07
Tags: drug response prediction, multi-omics integration, deep learning, translational, precision oncology

Proposes a structured multi-view approach (MulDR) to integrate multi-omics and drug structure for drug response prediction while preserving view-specific relational structure. This is relevant to translational modeling of therapy response/resistance, even if not pediatric-tumor-specific in the title/summary.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Drug response prediction (DRP) requires models that jointly characterize molecular drug structure and heterogeneous cellular states. Although recent deep learning approaches have integrated multi-omics profiles with drug representations, many treat omics views as parallel feature sources or merge them before relational modeling, which can weaken omics-view-specific cell-line relation structures. Here, we propose MulDR, a structur…

</details>

---

### [Emerging strategies for paediatric cancer predisposition: detection, surveillance and prevention](https://www.nature.com/articles/s41571-026-01217-z)
*Nature Reviews Clinical Oncology*  
Score: **0.69**
Published: 2026-11-06
Tags: pediatric oncology, cancer predisposition, surveillance, prevention, clinical

A clinical/translational perspective on identifying pediatric cancer predisposition syndromes to enable surveillance and prevention. While not computational, it is directly within pediatric oncology and may inform cohort stratification and germline-aware biomarker interpretation.

<details>
<summary>RSS summary</summary>

<p>Nature Reviews Clinical Oncology, Published online: 06 November 2026; <a href="https://www.nature.com/articles/s41571-026-01217-z">doi:10.1038/s41571-026-01217-z</a></p>Many childhood cancers occur in children with an underlying predisposition to developing cancer. Thus, identifying individuals with a cancer predisposition syndrome (CPS), ideally before the development of cancer, can enable regular clinical surveillance and potentially the introduction of preventive measures. In this Perspect…

</details>

---

### [Somatic structural variation discovery on cancer pangenome graphs](https://www.biorxiv.org/content/10.64898/2026.11.05.764632v1?rss=1)
*bioRxiv Genomics*  
Score: **0.66**
Published: 2026-11-06
Tags: structural variants, pangenome graphs, somatic SV calling, WGS, genomics method

Describes JANUS for calling somatic structural variants using a personal pangenome graph, aiming to reduce reference bias and computational costs versus linear-reference workflows. This is methodologically relevant for tumor WGS analyses where SVs are key drivers, though no pediatric context is stated here.

<details>
<summary>RSS summary</summary>

Somatic structural variants (SVs) drive tumorigenesis across most cancers, yet low variant allele frequencies in subclonal populations often necessitate multi-sample sequencing per patient. Existing somatic SV callers rely on whole-genome read alignments against a single linear reference, incurring high computational costs and reference bias. We present JANUS, which discovers somatic SVs directly within a personal pangenome, a colored de Bruijn-like graph of locally consistent parsing (LCP) core…

</details>

---

### [Detection and characterization of in cellulo protein–RNA complexity by irCLIP-RNP, irCLIPv2 and Re-CLIP](https://www.nature.com/articles/s41596-026-01459-2)
*Nature Protocols*  
Score: **0.63**
Published: 2026-11-06
Tags: protocol, CLIP, RNP, proteomics, transcriptomics

Protocol paper describing procedures to interrogate multi-protein ribonucleoprotein assemblies using transcriptomics and proteomics readouts. Useful for mechanistic biomarker programs involving RNA-binding proteins and post-transcriptional regulation, though not cancer-specific from the provided summary.

<details>
<summary>RSS summary</summary>

<p>Nature Protocols, Published online: 06 November 2026; <a href="https://www.nature.com/articles/s41596-026-01459-2">doi:10.1038/s41596-026-01459-2</a></p>Ribonucleoprotein assemblies are crucial regulators of RNA processing. Here, three interconnected procedures enabling the detection and interrogation of these multi-protein assemblies using transcriptomics and proteomics are presented.

</details>

---

### [Paxalisib with Trastuzumab for HER2-Positive Metastatic Breast Cancer with Active Brain Metastases: A Phase II Trial with Integrated Plasma and CSF ctDNA Analyses](https://aacrjournals.org/cancerrescommun/article/6/11/2619/789865/Paxalisib-with-Trastuzumab-for-HER2-Positive)
*AACR Cancer Research Communications*  
Score: **0.61**
Published: 2026-11-06
Tags: ctDNA, CSF, clinical trial, PI3K/mTOR, biomarkers, translational

A phase II trial combining a brain-penetrant PI3K/mTOR inhibitor with trastuzumab and integrating plasma/CSF ctDNA analyses, which is relevant to translational biomarker monitoring in the CNS setting. Cancer type is outside pediatric NB, but ctDNA-in-CSF concepts can transfer to pediatric neuro-oncology workflows.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle">Purpose:</div>Human epidermal growth factor receptor 2 (HER2)–positive breast cancer brain metastases (BCBM) frequently harbor PI3K/AKT/mTOR alterations. Preclinical HER2+ BCBM models demonstrated durable tumor regression with dual PI3K and mTOR inhibition. Our phase II trial evaluated paxalisib, a brain-penetrant dual PI3K/mTOR inhibitor, combined with trastuzumab.<div class="boxTitle">Patients and Methods:…

</details>

---

### [Targeting REV1 Exposes Context-Dependent Therapeutic Responses in Colorectal Cancer](https://www.biorxiv.org/content/10.64898/2026.11.04.764315v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.60**
Published: 2026-11-06
Tags: DNA damage tolerance, TLS, REV1, chemoresistance, organoids, preclinical

Studies REV1 inhibition in the DNA damage tolerance/TLS pathway and tests sensitivity in patient-derived colorectal cancer organoids. It is relevant to therapy resistance and preclinical modeling, though the specific disease context is CRC rather than pediatric cancers.

<details>
<summary>RSS summary</summary>

Translesion synthesis (TLS) promotes DNA damage tolerance and contributes to chemotherapy resistance by enabling replication across damaged DNA templates. REV1, a central TLS regulator, is a potential therapeutic target to enhance the efficacy of genotoxic chemotherapy. Here, we investigated the therapeutic consequences of pharmacologic and genetic REV1 targeting in colorectal cancer (CRC) and evaluated the sensitivity of patient-derived CRC organoids to the REV1 inhibitor JH-RE-06. JH-RE-06 enh…

</details>

---

### [Unveiling immune targets via TriStageHLA: high-throughput screening of peptide presentation, peptide-HLA binding, and immunogenicity](https://www.biorxiv.org/content/10.64898/2026.10.29.762831v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.55**
Published: 2026-11-09
Tags: neoantigen, HLA, immunogenicity, ESM-2, LoRA, immunotherapy

Proposes a three-stage epitope discovery pipeline including a LoRA-adapted ESM-2 classifier for early pruning plus HLA binding/presentation and immunogenicity prediction. Useful for immunotherapy target discovery, but not specifically tied to pediatric tumors or single-cell data in the provided summary.

<details>
<summary>RSS summary</summary>

High-throughput screening of HLA class I epitopes underpins neoantigen vaccine design, T cell therapies, and pathogen immune monitoring, yet proteome-scale prioritization remains challenging under severe class imbalance. We present TriStageHLA, a compact three-stage framework that decomposes epitope discovery into peptide-only screening (PRE), allele-specific binding/presentation prediction (BIND), and immunogenicity prediction (IM). PRE performs early pruning via a LoRA-adapted ESM-2 classifier…

</details>

---

### [Tumor metastasis co-opts liver hematopoiesis to promote immunosuppression](https://www.science.org/doi/abs/10.1126/science.aea3046?af=R)
*Science*  
Score: **0.52**
Published: 2026-11-05
Tags: metastasis, immunosuppression, tumor-immunology, mechanism

Addresses a tumor-driven mechanism of systemic immunosuppression via liver hematopoiesis in the context of metastasis. While not pediatric/NB-specific in the item text, immunosuppressive mechanisms can inform translational hypotheses and immune-biomarker interpretation.

<details>
<summary>RSS summary</summary>

Science, Volume 394, Issue 6824, November 2026. <br />

</details>

---

### [mFrag: Integrating External DNA Methylation References Improves Lung Cancer Early Detection by Fragmentomics](https://www.biorxiv.org/content/10.64898/2026.10.28.762499v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.44**
Published: 2026-11-09
Tags: liquid biopsy, cfDNA, fragmentomics, DNA methylation, early detection

Presents a liquid-biopsy method integrating methylation reference atlases with cfDNA fragmentomics for lung cancer detection. Methodologically relevant to noninvasive biomarker work, but the focus is early detection in lung cancer rather than pediatric oncology or treatment response.

<details>
<summary>RSS summary</summary>

Fragmentomics analyzes features of circulating cell-free DNA fragments (e.g., fragment length or end motifs) derived from whole-genome sequencing of liquid biopsy samples (blood or other body fluids), enabling non-invasive early cancer detection. Here we present mFrag, a framework that integrates external methylation reference atlases with fragmentomic features for lung cancer detection. mFrag comprises two complementary methods: mFrag-Weight, which uses DNA methylation information to weight fra…

</details>

---

### [Single-Cell Transcriptomics Reveals Tumorigenic Lineages and Interactions Persisting and Evolving From Premalignancy to Cancer in the Stomach](https://www.biorxiv.org/content/10.64898/2026.11.04.764362v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.43**
Published: 2026-11-06
Tags: scRNA-seq, tumor evolution, mouse model, gastric cancer, cell-cell interactions

Uses scRNA-seq (and high-parameter cytometry) to study lineage and interaction changes across premalignancy-to-cancer progression in a gastric cancer mouse model. It is adjacent single-cell cancer biology, but not directly tied to pediatric oncology or therapy response in the provided summary.

<details>
<summary>RSS summary</summary>

Background & Aims Gastric cancer is often diagnosed at advanced stages, making the premalignant window critical for early detection. However, the epithelial changes that characterize this stage remain poorly defined. We used CTLA4-knockdown (CTLA4KD) mice, which model human CTLA4 insufficiency-associated gastric cancer, that develops autoimmune gastritis, metaplasia, and adenocarcinoma, to define molecular changes across disease progression. Methods Single-cell RNA sequencing and 40-color spectr…

</details>

---

### [Genomic ascertainment refines the population prevalence and phenotype associations of germline TERT variants](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1012352)
*PLOS Genetics*  
Score: **0.38**
Published: 2026-11-05
Tags: germline, TERT, telomere biology, population genetics, predisposition

Focuses on population prevalence and phenotype associations of germline TERT variants, largely around telomere biology and non-pediatric outcomes in the summary. Potentially relevant background for germline predisposition considerations, but not central to translational pediatric tumor biology here.

<details>
<summary>RSS summary</summary>

<p>by Marena R. Niewisch, Sherifa A. Akinniyi, Kelvin C. de Andrade, Jung Kim, Logan P. Zeigler, Jeremy S. Haley, Gretchen Urban, David J. Carey, Douglas R. Stewart, Sharon A. Savage</p> Background <p>Individuals with deleterious germline variants in telomerase (encoded by <i>TERT</i>) are at high risk of pulmonary fibrosis (PF), liver disease (LD), bone marrow failure, head/neck squamous cell carcinoma (HNSCC), and other complications. Clinically ascertained cohorts are subject to biased ascert…

</details>

---

### [Transcriptome- and methylome-wide association studies prioritized lung cancer genes in East Asian populations](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1012332)
*PLOS Genetics*  
Score: **0.36**
Published: 2026-11-06
Tags: TWAS, MWAS, multi-omics, lung cancer, population genetics

Performs TWAS/MWAS to prioritize lung cancer genes with a focus on East Asian reference models and cross-population generalizability. It’s method-adjacent multi-omics association work, but not directly aligned with pediatric translational oncology or single-cell/resistance themes based on the summary.

<details>
<summary>RSS summary</summary>

<p>by Xingjie Hao, Minghan Qu, Yi Jiang, Cheng Wang, Ziwei Zhu, Pinpin Long, Chaolong Wang, Hongxia Ma, Shanshan Cheng</p> The relationship between gene expression, DNA methylation, and lung cancer remains elusive in East Asian (EAS) populations. Transcriptome- and methylome-wide association studies (TWAS/MWAS) help uncover gene–trait associations by integrating multi-omics data, yet models trained mainly on European references limit cross-population generalizability. To address this, we develop…

</details>

---

### [A Three-View System for CRISPR-Cas12a crRNA Activity Ranking: When to Ensemble, When Context Matters More](https://www.biorxiv.org/content/10.64898/2026.11.07.765023v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.34**
Published: 2026-11-09
Tags: CRISPR, Cas12a, guide-RNA-design, machine-learning, screening

Presents TriCasRank, a model combining sequence features, reference mapping, and historical activity summaries to rank Cas12a crRNA activity with cross-validation safeguards. Method is relevant to functional genomics screening design, but not specific to cancer, single-cell, or biomarkers in the provided summary.

<details>
<summary>RSS summary</summary>

CRISPR-Cas12a assays require effective crRNA selection, but experimentally screening all possible candidates is costly. We developed TriCasRank, a ranking model that combines sequence features with EasyDesign reference mapping and activity summaries from related training records. Historical summaries were calculated within each training fold to avoid information leakage. In five-fold target-grouped cross-validation using out-of-fold predictions, TriCasRank improved the Spearman correlation from …

</details>

---

### [Integrated genotyping, annotation and epigenetic characterization of transposable elements with GraffiTE](https://www.biorxiv.org/content/10.64898/2026.11.06.764993v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.31**
Published: 2026-11-09
Tags: long-read, methylation, transposable-elements, Nextflow, pipeline

Introduces a Nextflow pipeline to genotype structural variants/mobile element insertions and retain allele-specific CpG methylation information from long-read sequencing. Strongly computational, but TE-focused rather than directly aligned to cancer biomarker/therapy response in the item text.

<details>
<summary>RSS summary</summary>

Genotyping, annotating and characterizing TE variation leaped forward with long-read sequencing (LRS) and de novo assembly. LRS techniques also reveal CpG methylation, which enables measuring the impact of TEs on epigenomes and sheds light on the role of TEs as regulatory sequences. Yet, TE genotyping tools discard the valuable allele-specific epigenetic information. We introduce GraffiTE v1.1, the first Nextflow pipeline that genotypes SVs, annotates polymorphic mobile elements (MEI), and maps …

</details>

---

### [Membrane tension and cell–nucleus mechanical coupling regulate confined cell migration](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014909)
*PLOS Comp Bio*  
Score: **0.24**
Published: 2026-11-05
Tags: cell-migration, mechanobiology, computational-modeling

Computational/quantitative focus on mechanics of confined cell migration and cell–nucleus coupling. Could be tangentially relevant to metastasis biology, but no biomarker, omics, or pediatric oncology connection is indicated in the summary snippet.

<details>
<summary>RSS summary</summary>

<p>by David Hernandez-Aristizabal, Anotida Madzvamuse, Rachele Allena</p> Cell migration is a fundamental physiological process essential from embryogenesis to cancer metastasis. It is primarily governed by the coupled mechanics of the plasma membrane (PM) and the nucleus, together with cytoskeletal dynamics. Being the largest and stiffest organelle, the nucleus is considered a rate-limiting factor in confined migration. However, nuclear positioning varies across cell types and migration modes, …

</details>

---

### [BEHAV3D-Explorer enables natural language-assisted profiling of 3D cellular dynamics](https://www.biorxiv.org/content/10.64898/2026.11.06.764820v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.22**
Published: 2026-11-09
Tags: 3D microscopy, image analysis, natural language interface, cell dynamics, tooling

Describes an interactive platform for analyzing live 3D microscopy with natural-language assistance for workflow configuration. Useful tooling, but it is not clearly tied to omics, biomarker discovery, or pediatric oncology in the provided summary.

<details>
<summary>RSS summary</summary>

Obtaining quantitative measurements from live 3D microscopy remains difficult to implement across research settings, particularly without profound image-analysis expertise. We present BEHAV3D-Explorer, an interactive, user-friendly platform that uses natural-language guidance to help users select analysis strategies, tune parameters and configure workflows. BEHAV3D-Explorer enables time-resolved behavioral classification and integrates it with morphological, interaction and functional analyses, …

</details>

---

### [Air Pollution and Atopic Dermatitis in Children](https://jamanetwork.com/journals/jamapediatrics/fullarticle/2856356)
*JAMA Pediatrics*  
Score: **0.03**
Published: 2026-11-09
Tags: pediatrics, epidemiology, meta-analysis, dermatology

Systematic review/meta-analysis on air pollutants and atopic dermatitis prevalence in children. Pediatric but not related to oncology, omics, or translational cancer biomarker work.

<details>
<summary>RSS summary</summary>

This systematic review and meta-analysis examines the associations between ambient air pollutants and atopic dermatitis prevalence in children.

</details>

---

### [Satisfaction, Regret, and Harm After Hormone Therapy in Youth](https://jamanetwork.com/journals/jamapediatrics/fullarticle/2856355)
*JAMA Pediatrics*  
Score: **0.02**
Published: 2026-11-09
Tags: pediatrics, cohort-study, health-outcomes

Cohort study on outcomes after hormone therapy initiated in youth. Not related to pediatric oncology, computational biology for omics, or translational cancer research.

<details>
<summary>RSS summary</summary>

This cohort study examines satisfaction, regret, and perceived harm among individuals who initiated hormone therapy as youth and identifies demographic and clinical characteristics associated with regret or harm and current hormone therapy status.

</details>

---

### [Antibacterial properties of commercially available plant extracts containing tannins against Escherichia coli](https://www.nature.com/articles/s41598-026-79871-7)
*Scientific Reports*  
Score: **0.00**
Published: 2026-11-09
Tags: antibacterial, microbiology, natural-products

Tests antibacterial activity of tannin-containing plant extracts against E. coli. Not relevant to pediatric translational oncology or computational omics.

<details>
<summary>RSS summary</summary>

<p>Scientific Reports, Published online: 09 November 2026; <a href="https://www.nature.com/articles/s41598-026-79871-7">doi:10.1038/s41598-026-79871-7</a></p>Antibacterial properties of commercially available plant extracts containing tannins against <i>Escherichia coli</i>

</details>

---

## AI (2 shown / 2 total)

### [BioMeme: biological response memory for virtual cells](https://www.biorxiv.org/content/10.64898/2026.11.07.765043v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.71**
Published: 2026-11-09
Tags: perturbation modeling, drug response, interpretability, virtual cells, representation learning

Introduces a “biological response memory” approach for virtual cells that predicts perturbation responses while recalling supporting experiments rather than relying purely on absorbed neural network weights. This is relevant to interpretable perturbation/response modeling workflows for drug mechanism and combination prioritization.

<details>
<summary>RSS summary</summary>

Virtual cells predict responses to unmeasured perturbations, but neural models absorb measurements into weights, require retraining and cannot cite supporting experiments. Here we introduce biological response memory, replacing a learned state-to-response layer with a memory of measured experiments. BioMeme predicts next-dose responses in a new cell type from measured lower-dose states, recalling similarly responding compounds in that cell and the same compound in a reference cell line. Each pre…

</details>

---

### [TriTower-m6Am: A triple-tower heterogeneous deep learning architecture integrating semantic, sequential, and structural information for mRNA N6,2′-O-dimethyladenosine site prediction](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014939)
*PLOS Comp Bio*  
Score: **0.22**
Published: 2026-11-05
Tags: deep-learning, RNA-modifications, m6Am, prediction

Deep learning model for predicting m6Am RNA modification sites by integrating multiple sequence/structure representations. This is primarily a general RNA-modification prediction method with no explicit link to cancer, single-cell, or translational biomarker applications in the provided text.

<details>
<summary>RSS summary</summary>

<p>by Kaifeng Xiong, Jianhua Jia</p> Background <p>N<sup>6</sup>,2′-O-dimethyladenosine (m6Am) is a cap-proximal mRNA modification deposited by PCIF1 (phosphorylated C-terminal domain (CTD)-interacting factor 1) at the first transcribed nucleotide of eukaryotic mRNAs. Knowing where m6Am sites sit would help explain how cells tune mRNA stability and translation, but current predictors typically depend on a single representation and do not jointly model multiple views of the same RNA sequence.</p>…

</details>

---

## Methods (3 shown / 3 total)

### [spaDWGCN: a spatial distance weighted graph convolution network for accurate spatial tissue structure identification](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag623/9022190?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.58**
Published: 2026-11-08
Tags: spatial transcriptomics, graph neural networks, tissue architecture, GCN, microenvironment

Introduces a spatial distance-weighted GCN to identify spatial tissue architecture using expression plus coordinates, addressing limitations of binary adjacency. This is applicable to tumor microenvironment niche detection in spatial datasets, though the summary does not indicate clinical/translational evaluation.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>The continuous development of spatial transcriptomics enables us to obtain gene expression matrices alongside spatial coordinate labels for cells or spots, providing tools for understanding spatial tissue microenvironments and spatially resolved cell–cell interactions. Although numerous methods have been developed to identify spatial tissue architecture, most of them adopt binary adjacency matrices that treat all neighboring cell…

</details>

---

### [MultiMHC: a multiscale learning framework integrating sequence and residue contact topology for peptide–MHC binding prediction](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag466/9017776?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.46**
Published: 2026-11-06
Tags: peptide-MHC, binding prediction, deep learning, immunotherapy, epitopes

Develops a deep learning framework that incorporates sequence and residue contact topology to improve peptide–MHC binding prediction. Relevant for cancer immunotherapy pipelines, but peripheral to the user’s core single-cell/multi-omics resistance/biomarker focus.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Accurate prediction of peptide–major histocompatibility complex (pMHC) binding is central to understanding antigen presentation and identifying T-cell epitopes in cancer immunotherapy. Although recent methods have achieved notable progress, most approaches rely predominantly on sequence-based features and only partially capture structural determinants of pMHC recognition. In this study, we present MultiMHC, a multiscale deep lear…

</details>

---

### [Correction: Deep learning-based prediction of lymph node metastasis and occult tumor cells in gastric cancer using histopathological images: a retrospective study](https://www.nature.com/articles/s41416-026-03640-6)
*British Journal of Cancer*  
Score: **0.00**
Published: 2026-11-09
Tags: correction, digital-pathology, deep-learning, gastric-cancer

This item is a correction notice for a gastric cancer histopathology deep learning study, not new methods/results in the RSS entry. Also not aligned with pediatric oncology or single-cell/multi-omics workflows based on the provided text.

<details>
<summary>RSS summary</summary>

<p>British Journal of Cancer, Published online: 09 November 2026; <a href="https://www.nature.com/articles/s41416-026-03640-6">doi:10.1038/s41416-026-03640-6</a></p>Correction: Deep learning-based prediction of lymph node metastasis and occult tumor cells in gastric cancer using histopathological images: a retrospective study

</details>

---

## Other (7 shown / 7 total)

### [scYOU: a multi-scale yoked optimization and unsupervised framework for analysis of cellular heterogeneity and complex-associated protein module recovery](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag626/9022372?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.82**
Published: 2026-11-09
Tags: proteomics, cellular heterogeneity, protein modules, unsupervised learning, multi-scale, computational method

Presents an unsupervised, multi-scale framework intended to jointly model cell-level heterogeneity and protein-level module structure in high-resolution proteomics. This is directly relevant to tumor state discovery and protein program inference from (single-cell or high-resolution) proteomic data.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Advances in high-resolution proteomic technologies are enabling systematic characterization of cellular heterogeneity and protein functional programs at increasing depth. However, existing computational approaches typically focus on either cell-level organization or protein-level structure, lacking a unified framework that jointly models both entities. Here we present scYOU, a multi-scale yoked optimization and unsupervised frame…

</details>

---

### [Contrasting within-donor gene covariation identifies disease-associated cells and genes in case-control single-cell studies](https://www.biorxiv.org/content/10.64898/2026.11.06.764894v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.80**
Published: 2026-11-09
Tags: scRNA-seq, case-control, donor effects, disease-associated cells, gene covariation, computational method

Proposes IDC to detect disease-associated subpopulations in case-control scRNA-seq by leveraging within-donor heterogeneity and coordinated gene covariation, addressing dilution of signals and donor/batch confounding. This directly matches needs in clinical single-cell studies where only subsets of cells are affected (e.g., resistant states).

<details>
<summary>RSS summary</summary>

In case-control single-cell studies, disease often affects only a subset of cells of a given type within each patient. Comparing all cells from patients with all cells from control donors therefore dilutes disease signals, while donor variation and recorded or unrecorded batch effects can be confounded with disease status. We present IDC (Identification of Disease-associated Cells), which uses this cellular heterogeneity as signal. A mixture of affected and unaffected cells creates coordinated v…

</details>

---

### [Accelerating discovery and integration from RNA to Protein with R2T2P](https://www.biorxiv.org/content/10.64898/2026.11.05.764226v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.78**
Published: 2026-11-09
Tags: transcriptomics, proteomics, RNA-to-protein, integration

Integrates transcriptomics with proteomics translation pipeline, enabling unified RNA-to-protein analyses with clear translational potential.

---

### [SPACE-Tag Enables High-Throughput Spatial Chromatin Profiling With CUT&Tag](https://www.biorxiv.org/content/10.64898/2026.11.06.764689v1?rss=1)
*bioRxiv Genomics*  
Score: **0.72**
Published: 2026-11-07
Tags: spatial-omics, CUT&Tag, epigenomics, multi-omics, methods

Presents a spatial epigenomics CUT&Tag framework designed to integrate epigenetic marks with standard spatial transcriptomics workflows, supporting multi-omic mapping of regulatory landscapes in tissue. The focus on scalable spatial chromatin profiling is methodologically transferable to tumor microenvironment and cell-state studies.

<details>
<summary>RSS summary</summary>

Spatial epigenomics maps chromatin states within intact tissues, yet current methods remain limited by accessibility and scale. We developed SPACE-Tag, a spatial CUT&amp;Tag framework that converts epigenetic marks-associated DNA into polyadenylated RNA through T7 based linear amplification, enabling direct integration with standard spatial transcriptomic workflows. Applied to the mouse brain, SPACETag resolved activating and repressive chromatin landscapes, recovered region-specific gene regula…

</details>

---

### [scHiR learns joint cellular representations from paired single-cell transcriptome and 3D genome data](https://www.biorxiv.org/content/10.64898/2026.11.06.764816v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.70**
Published: 2026-11-09
Tags: single-cell multi-omics, scRNA-seq, scHi-C, representation learning, 3D genome

Focuses on integrative representation learning for matched scRNA-seq and scHi-C measurements within the same cells to capture complementary state information. This is a broadly transferable multi-omics integration method for regulatory-state inference, though not explicitly cancer-focused in the provided text.

<details>
<summary>RSS summary</summary>

Simultaneously sequenced scRNA-seq and scHi-C technologies provide matched measurements of gene expression and three-dimensional (3D) genome organization within individual cells, offering complementary insights into cellular state. However, integrative analysis remains challenging due to the heterogeneous structures and characteristics of these modalities. To realize the full potential of such data, we need effective computational methods to learn joint cellular representations. Here, we present…

</details>

---

### [CalDyn: an advanced high-throughput automated algorithm for monitoring intracellular calcium dynamics and electrophysiological signals in patient-derived cells](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag602/9022079?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.43**
Published: 2026-11-07
Tags: software, calcium-imaging, electrophysiology, patient-derived, QC

Introduces a MATLAB-based tool for automated peak detection, waveform characterization, event classification, and QC for calcium dynamics/electrophysiology data in patient-derived cells. Potentially useful if your preclinical models include functional phenotyping, but not directly tied to omics or oncology in the summary.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle">Background</div>Quantitative analysis of intracellular calcium and electrophysiological signals is challenging when signals exhibit baseline drift, low-amplitude events, irregular waveforms, or multicomponent transients. We developed CalDyn, a MATLAB-based graphical application for automated peak detection, waveform characterization, event classification, batch processing, quality control and data export.<di…

</details>

---

### [Evidence and measures of the digital usage divide from 40 million Windows devices](https://www.nature.com/articles/s41467-026-78689-1)
*Nature Communications*  
Score: **0.01**
Published: 2026-11-09
Tags: digital-divide, population-data

Develops measures of digital engagement using telemetry from Windows devices to characterize digital inequality. Outside biomedical and oncology scope.

<details>
<summary>RSS summary</summary>

<p>Nature Communications, Published online: 09 November 2026; <a href="https://www.nature.com/articles/s41467-026-78689-1">doi:10.1038/s41467-026-78689-1</a></p>Using usage data from 40 million Windows devices, this paper develops measures of digital engagement in the USA, revealing digital inequalities that are not captured by traditional measures of internet access alone.

</details>

---
