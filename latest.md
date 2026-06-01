# Weekly ToC Digest (week of 2026-06-01)

This digest is automatically generated from this week's RSS items and categorized into Neuroblastoma, AI, Methods, and Other.

**Included:** 33 (all ranked papers)  
**Scored:** 33 total items

**Models:** `gpt-5-nano` first-pass, `gpt-5.2` re-rank top 40

---

## Neuroblastoma (19 shown / 19 total)

### [A Foundation Model for the Cancer Genome](https://www.biorxiv.org/content/10.64898/2026.05.27.728319v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.83**
Published: 2026-06-01
Tags: foundation model, cancer genomics, self-supervised learning, copy-number, biomarkers, translational

Proposes a self-supervised foundation model trained on large-scale tumor genomics to interpret co-occurring alterations and genome-wide copy-number patterns, directly aligning with scalable biomarker/response modeling in cancer. Strongly transferable to translational oncology workflows even without pediatric/neuroblastoma specificity in the title/summary.

<details>
<summary>RSS summary</summary>

Cancer is a disease of the genome, in which somatic mutations and copy-number alterations determine tumour identity, clinical behaviour, and response to therapy. Consortium-scale sequencing has profiled hundreds of thousands of tumours, yet clinical interpretation still proceeds one alteration at a time against hand-curated knowledgebases, often ignoring co-occurring alterations and the genome-wide copy-number pattern. Self-supervised foundation models pretrained on unlabelled corpora have produ…

</details>

---

### [Genomic-Adjusted Radiation Dose from Bulk RNA Sequencing for Personalized Radiotherapy](https://www.biorxiv.org/content/10.64898/2026.05.29.728725v1?rss=1)
*bioRxiv Genomics*  
Score: **0.74**
Published: 2026-05-30
Tags: radiotherapy, bulk RNA-seq, treatment response, predictive biomarker, clinical stratification, translational

Links bulk RNA-seq to a radiosensitivity-based metric (GARD/RSI) intended to personalize radiotherapy dose and predict benefit, fitting treatment-response biomarker/stratification themes. While not pediatric-specific, the transcriptomics-to-clinical-decision framing is highly translational.

<details>
<summary>RSS summary</summary>

Radiotherapy is delivered to more than half of all patients with cancer yet is prescribed using uniform physical doses despite well-established interpatient variability in biological response. The genomic-adjusted radiation dose (GARD), derived from the radiosensitivity index (RSI), integrates tumor transcriptomics with radiation dose to estimate patient-specific treatment effect, and has been clinically validated as a predictor of radiotherapy benefit across diverse disease sites, including bre…

</details>

---

### [TF-DWGNet: a directed weighted graph neural network with tensor fusion for multi-omics cancer subtype classification](https://academic.oup.com/nargab/article/doi/10.1093/nargab/lqag054/8698434?rss=1)
*NAR Genomics & Bioinformatics*  
Score: **0.71**
Published: 2026-05-30
Tags: multi-omics, graph neural network, cancer subtypes, data integration, machine learning

Introduces a GNN + tensor fusion framework to integrate heterogeneous multi-omics for cancer subtype classification, matching multi-omics integration and robust computational modeling priorities. Applicability is broad and could transfer to pediatric tumor stratification settings.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Integration and analysis of multi-omics data provide valuable insights for improving cancer subtype classification. However, such data are inherently heterogeneous, high-dimensional, and exhibit complex intra- and inter-modality dependencies. Graph neural networks provide a principled framework for modeling these structures, but existing approaches often rely on prior knowledge or predefined similarity networks that produce eithe…

</details>

---

### [Globular domain histone H3R131C mutation remodels chromatin accessibility to promote oncogenic transcriptional programs](https://www.biorxiv.org/content/10.64898/2026.05.27.728333v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.68**
Published: 2026-05-29
Tags: chromatin, epigenetics, oncogene

Epigenetic remodeling by an oncohistone offers mechanistic insight into transcriptional programs; has translational relevance for epigenetic therapy and biomarker discovery.

---

### [Personalized reference genome-based pipeline reveals comprehensive haplotype-resolved views of cancer genomes](https://www.biorxiv.org/content/10.64898/2026.05.28.728591v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.67**
Published: 2026-05-30
Tags: cancer genomics, haplotype-resolved, structural variants, copy-number, DNA methylation, pipeline

Presents a pipeline to call somatic SNVs/SVs/CNAs and methylation using personalized diploid references and haplotype resolution, addressing hard-to-detect alterations relevant to actionable genomics. Method is directly transferable to clinical/translational cancer genome analysis pipelines.

<details>
<summary>RSS summary</summary>

Cancer genome analysis relies on standard human reference genomes but detecting somatic alterations in highly repetitive or individual-specific regions remains challenging. We developed the Personalized Reference genome-based Cancer Genome Analysis Pipeline (PRCGAP), to our knowledge, the first comprehensive pipeline integrating haplotype-resolved analyses of somatic point mutations, structural variants, copy number, and DNA methylation on personalized diploid reference genomes. We applied PRCGA…

</details>

---

### [Rare RNA Polymerase II failure modes mark the cancer-driving genes most affected by epigenetic perturbation](https://www.biorxiv.org/content/10.64898/2026.05.28.728581v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.66**
Published: 2026-05-31
Tags: transcription, RNA Pol II, epigenetics, perturbation, computational method, ChIP-seq

Introduces a computational tool to compare Pol2 activity state shifts using ChIP data before/after chemical perturbation, expanding beyond promoter-proximal pausing. The focus on perturbation response and identifying cancer-driving genes affected by epigenetic perturbation is methodologically relevant even outside pediatrics.

<details>
<summary>RSS summary</summary>

RNA Polymerase II (Pol2) transcribes genes through a complex life cycle (initiation, pausing, elongation, co-transcriptional splicing, termination, and recycling). Chromatin immunoprecipitation of Pol2 before and after chemical perturbation has identified promoter-proximal accumulation (pausing) as a critical step in the transcription genome-wide. However, the full landscape of Pol2 responses has not been well characterized. Here, we introduce a tool for comparing Pol2 Activity State Shifts (com…

</details>

---

### [Assessing and Optimizing Low-Frequency Somatic Mutation Detection: A Multi-Platform High-Throughput Sequencing Perspective](https://www.biorxiv.org/content/10.64898/2026.05.28.728367v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.63**
Published: 2026-06-01
Tags: low VAF, somatic variants, cfDNA, benchmarking, sequencing platforms, MRD

Systematically compares sequencing platforms for low-VAF somatic variant detection using GIAB and cfDNA reference standards, which is relevant to ctDNA/MRD-style assay design and benchmarking. Emphasis is methodological/technical but directly impacts translational biomarker detection sensitivity.

<details>
<summary>RSS summary</summary>

The availability of multiple commercial short-read sequencing platforms necessitates systematic cross-platform performance comparisons, particularly for challenging applications such as low-frequency somatic mutation detection. Here, a large-scale targeted sequencing dataset from five Genome in a Bottle (GIAB) human genomic DNA reference standards, HG001 to HG005, alongside Twist Biosciences cfDNA reference standards featuring 1% variant allele frequency (VAF), was generated by six platforms (No…

</details>

---

### [Auranofin potentiates cisplatin response through context-dependent NOTCH-associated signaling states in endometrial cancer](https://www.biorxiv.org/content/10.64898/2026.05.27.728338v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.63**
Published: 2026-05-31
Tags: therapy response, drug resistance, cisplatin, NOTCH, combination therapy, translational

Directly addresses therapeutic resistance and drug response modulation (auranofin altering platinum responsiveness) tied to NOTCH-associated signaling states. While not pediatric or neuroblastoma, the response/resistance framing and signaling-state context are translationally aligned.

<details>
<summary>RSS summary</summary>

Therapeutic resistance remains a major challenge in advanced and recurrent endometrial cancer (EC). Aberrant NOTCH signaling has been associated with aggressive tumor behavior and therapeutic resistance across multiple malignancies, yet its therapeutic significance in EC remains incompletely defined. We investigated whether auranofin (AuR), a noncanonical modulator of NOTCH signaling through the transcriptional effector RBPJ, alters platinum responsiveness in EC models. Elevated NOTCH3 copy-numb…

</details>

---

### [Multiomic dissection of HR+/HER2- invasive lobular breast carcinoma reveals mobilized yet dysfunctional anti-tumor immunity shaped by tumor-stroma crosstalk and impaired antigen presentation](https://www.biorxiv.org/content/10.64898/2026.05.28.728418v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.60**
Published: 2026-05-29
Tags: multi-omics, tumor microenvironment, immunology, antigen presentation, tumor-stroma crosstalk, translational

Uses a multi-omics approach to characterize immune dysfunction, tumor–stroma crosstalk, and impaired antigen presentation in a clinically defined subtype, which is conceptually transferable to tumor microenvironment analyses. The clinical-translation angle (ICI context) supports relevance despite cancer-type mismatch.

<details>
<summary>RSS summary</summary>

Context: Immunotherapy based on immune checkpoint inhibitors (ICI) revolutionized the treatment of triple-negative (TN) breast carcinomas (BC), but remains more challenging in HR+/HER2- BCs. Because invasive lobular carcinomas (ILC) generally exhibit low immune infiltration, ICIs were largely overlooked in this pathological type. The only clinical trial of ICIs dedicated to ILCs showed disappointing results, notably in HR+/HER2- cases. The immune landscape of HR+/HER2- ILCs has been poorly descr…

</details>

---

### [Ultra-efficient High Resolution 3D Reconstruction of Spatial Omics Data with Neural Transcriptomic Field](https://www.biorxiv.org/content/10.64898/2026.05.28.726140v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.57**
Published: 2026-06-01
Tags: spatial omics, 3D reconstruction, neural field, representation learning, tissue architecture

Addresses 3D reconstruction of spatial omics volumes from sparse 2D sections using a neural field approach, potentially useful for tumor architecture and niche analysis. Relevance is methodological (spatial) rather than directly translational/clinical in the provided summary.

<details>
<summary>RSS summary</summary>

Biological tissues are inherently three-dimensional (3D) ecosystems where spatial architecture dictates cellular function. While spatial omics technologies have revolutionized molecular profiling, they are largely restricted to isolated two-dimensional (2D) tissue sections. Existing computational methods attempting to reconstruct 3D volumes from sparse slices rely heavily on local slice-to-slice interpolation, struggling to balance high-fidelity reconstruction, noise reduction, and atlas-scale e…

</details>

---

### [Integrative prioritization of clinically and biologically relevant long noncoding RNAs across gastrointestinal cancers](https://www.biorxiv.org/content/10.64898/2026.05.26.728026v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.52**
Published: 2026-05-29
Tags: lncRNA, biomarkers, integration, bulk RNA-seq, single-cell, pan-cancer

Develops an integrative framework using bulk and single-cell transcriptomic resources to classify/prioritize lncRNAs across multiple cancers. The integrative prioritization concept could transfer to biomarker discovery pipelines, though the application is GI-specific.

<details>
<summary>RSS summary</summary>

Across gastrointestinal (GI) cancers, shared malignant programs are layered onto strong anatomical, lineage, and microenvironmental variation, making it difficult to distinguish disease-relevant long noncoding RNAs (lncRNAs) from context-dependent transcriptional signals. We developed a pan-GI integrative framework to classify lncRNAs across colorectal adenocarcinoma, gastric adenocarcinoma, and esophageal cancer using bulk and single-cell transcriptomic resources. This framework evaluates lncRN…

</details>

---

### [TNBC Spatial Transcriptomic Analysis across Clinical States Reveals Subtype-Specific Networks and Immunosuppressive Niches](https://aacrjournals.org/cancerrescommun/article/6/5/1246/785468/TNBC-Spatial-Transcriptomic-Analysis-across)
*AACR Cancer Research Communications*  
Score: **0.39**
Published: 2026-05-29
Tags: spatial transcriptomics, tumor microenvironment, networks, immunosuppression, TNBC

Uses spatial transcriptomics across clinical states in TNBC to derive subtype-specific networks and immunosuppressive niches, which is relevant conceptually to resistance/microenvironment. However, it is breast-cancer specific and the method details/transferability are limited in the provided summary.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle"></div>Breast cancer is a heterogeneous disease composed of distinct molecular subtypes that influence prognosis and treatment response, with subtype discordance between primary and metastatic tumors contributing to therapeutic failure. Using GeoMx spatial transcriptomics, we profiled nonmetastatic primary tumors, metastatic primary tumors, and lymph node (LN) metastases to characterize transcriptional and im…

</details>

---

### [Targeting CD73-A2aR-Mediated Adenosine Signaling at the Tumor-Immune Interface Overcomes Radioresistance](https://www.biorxiv.org/content/10.64898/2026.05.26.727904v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.33**
Published: 2026-05-29
Tags: radioresistance, CD73, adenosine signaling, tumor-immune interface, breast cancer

Focuses on mechanistic immuno-metabolic checkpoint biology (CD73/adenosine/A2aR) affecting radioresistance in breast cancer models. Potentially relevant to resistance biology, but appears primarily wet-lab/therapy-mechanism rather than computational omics.

<details>
<summary>RSS summary</summary>

Background: Radiotherapy efficacy is constrained by an immunosuppressive tumor microenvironment (TME) enriched in extracellular adenosine and suppressive myeloid populations that attenuate cytotoxic T-cell responses. The CD73-adenosine-A2a/A2b receptor axis represents a key metabolic immune checkpoint; however, the relative contributions of tumor cell intrinsic versus host-derived adenosine signaling to radiotherapy response remain incompletely defined. Methods: Using orthotopic murine breast ca…

</details>

---

### [Cutaneous inflammation accelerates the premalignant expansion of melanocytes bearing oncogenic mutations](https://www.biorxiv.org/content/10.64898/2026.05.27.728115v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.24**
Published: 2026-05-29
Tags: melanoma, inflammation, microenvironment, mouse model, premalignant

Focuses on immune perturbations in a mouse model driving premalignant melanocyte expansion, emphasizing microenvironmental inflammation. The computational/omics angle is not evident from the provided summary, and it is not aligned with pediatric translational oncology priorities.

<details>
<summary>RSS summary</summary>

How the cutaneous microenvironment influences early melanomagenesis is poorly understood. Here, we assessed the effects of three immune perturbations on premalignant melanocyte expansion in an autochthonous mouse model of disease. Depletion of regulatory T (Treg) cells markedly accelerated melanoproliferation, an unexpected phenotype that was associated with monocyte and macrophage infiltration, the production of inflammatory and angiogenic factors, and vascular leakage. In line with these obser…

</details>

---

### [Conformational drugging pockets at the surface of c-MYC proto-oncongenes](https://www.biorxiv.org/content/10.64898/2026.05.28.728384v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.23**
Published: 2026-05-29
Tags: MYC, druggability, structural biology, oncogene

Discusses structural/druggability aspects of MYC, an important oncogene, but the snippet reads as structural biology/drug discovery rather than omics, single-cell, or translational biomarker work. Relevance is indirect to computational oncology priorities.

<details>
<summary>RSS summary</summary>

The MYC proto-oncogene encodes a master transcription factor essential for regulating cellular homeostasis, growth, and metabolism. However, its deregulation is involved in up to 70% of human cancers, where it drives uncontrolled cellular proliferation. Structurally, MYC is characterized by an N-terminal transactivation domain and a C-terminal bHLHLZ domain critical for DNA binding and dimerization with MAX. Despite its clinical significance, MYC has remained historically "undruggable" due to it…

</details>

---

### [Unveiling the novel role of PGAM5 in rewiring metabolism through PI3K/AKT/mTOR signaling in acute myelogenous leukemia](https://www.nature.com/articles/s41598-026-55582-x)
*Scientific Reports*  
Score: **0.18**
Published: 2026-06-09
Tags: AML, metabolism, PI3K-AKT-mTOR, mechanism

Acute myelogenous leukemia metabolism/signaling study; computational/translational biomarker components are not indicated in the short summary line. Outside core neuroblastoma/omics-method focus based on available information.

<details>
<summary>RSS summary</summary>

<p>Scientific Reports, Published online: 09 June 2026; <a href="https://www.nature.com/articles/s41598-026-55582-x">doi:10.1038/s41598-026-55582-x</a></p>Unveiling the novel role of PGAM5 in rewiring metabolism through PI3K/AKT/mTOR signaling in acute myelogenous leukemia

</details>

---

### [Population-specific heterogeneity in ontogeny of the broadly-conserved blood transcriptional program during the first week of life](https://www.nature.com/articles/s41467-026-73244-4)
*Nature Communications*  
Score: **0.15**
Published: 2026-06-01
Tags: neonatal, immunology, transcriptomics, population differences

Studies neonatal immune development across populations, relevant to pediatrics but not oncology, biomarkers for cancer, or therapy response. Any transfer to tumor immunology is indirect based on the provided abstract.

<details>
<summary>RSS summary</summary>

<p>Nature Communications, Published online: 01 June 2026; <a href="https://www.nature.com/articles/s41467-026-73244-4">doi:10.1038/s41467-026-73244-4</a></p>This study compares immune development (ie. ontogeny) in healthy newborns across two distinct populations, revealing differences on top of a strong core pattern of development. These differences can inform the influence of genetic or environmental factors on newborn immunity

</details>

---

### [Fully synthetic replication of complex real biological cell clusters using a novel cluster-based ‘Rosetta-Routine’ computational modelling process](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014280)
*PLOS Comp Bio*  
Score: **0.12**
Published: 2026-05-29
Tags: flow cytometry, synthetic data, computational modeling, cluster simulation

Focuses on synthetic replication of flow cytometry cell clusters for modeling/diagnostics support. It is computational but not clearly tied to single-cell transcriptomics, multi-omics integration, or oncology biomarker discovery in the provided summary.

<details>
<summary>RSS summary</summary>

<p>by Bradley Mason, Laura Justham, Liam Whitby, Alison Whitby, Stuart Scott, Samuel Nti, Jon Petzing</p> Flow cytometry (FC) is essential for the precise quantification and characterisation of individual cell populations in a larger heterogenous cell suspension. FC analysis provides a foundation for advanced clinical diagnostics and is a key component in many life-saving therapeutic strategies across a broad range of medical conditions. However, clinical, industrial and research laboratories al…

</details>

---

### [Graph-based pangenome provides insights into the adaptive evolution of Cucurbita pepo](https://www.biorxiv.org/content/10.64898/2026.05.27.728270v1?rss=1)
*bioRxiv Genomics*  
Score: **0.02**
Published: 2026-05-29
Tags: pangenome, plant genomics, graph genome, evolution

Plant pangenome/evolution study (Cucurbita pepo) is outside pediatric oncology and computational single-cell/biomarker priorities. Limited methodological transfer is evident from the snippet.

<details>
<summary>RSS summary</summary>

Understanding how crops respond to environmental variation is crucial for biodiversity conservation and food security. Cucurbita pepo (pumpkin, squash, gourd) is one of the first domesticated crop species and exhibits remarkable phenotypic and ecological diversity, making it a powerful system for investigating the genomic basis of adaptation. Here, we constructed a graph-based C. pepo pangenome using nine chromosome-level assemblies and identified 229,431 high-confidence structural variants (SVs…

</details>

---

## AI (0 shown / 0 total)

## Methods (5 shown / 5 total)

### [TRACE: a graph-based workflow for TCR-epitope prioritization and tumor-reactive T-cell identification](https://www.biorxiv.org/content/10.64898/2026.05.27.728217v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.72**
Published: 2026-05-31
Tags: immuno-oncology, TCR, epitope, graph methods, single-cell, tumor-reactive T cells

Presents a graph-based computational workflow to prioritize TCR–epitope interactions and identify tumor-reactive T cells, a common need in single-cell immuno-oncology analyses. The method focus (graph modeling over independent pairs) appears broadly transferable to scTCR-seq + transcriptomics workflows.

<details>
<summary>RSS summary</summary>

Accurate prioritization of T-cell receptor (TCR)-epitope interactions and identification of tumor-reactive T cells are important but difficult steps in immunotherapy-oriented bioinformatics workflows. Existing methods typically address these tasks separately and either model TCR-epitope pairs as independent observations or rely primarily on transcriptomic signatures. In this study, we present TRACE (TCR-epitope pRioritization And T-Cell idEntification), a graph-based computational workflow that …

</details>

---

### [Genotype and methylation interact to reconfigure transcriptional regulation in colorectal cancer](https://www.biorxiv.org/content/10.64898/2026.05.27.728350v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.50**
Published: 2026-05-30
Tags: regulatory genomics, DNA methylation, eQTL-like, multi-omics, colorectal cancer

Analyzes how methylation context modifies genotype–expression regulatory relationships using paired tumor/NAT samples, aligning with regulatory inference and multi-omics (genotype/epigenome/expression). Cancer-type specific (CRC), but concepts may transfer to pediatric tumor regulatory studies.

<details>
<summary>RSS summary</summary>

Background Transcriptional regulation is shaped by both genomic variants and the environment. Yet, how the regulatory effects of genomic variants are reconfigured by dynamic epigenomic changes during tumorigenesis remains incompletely understood. Methods We investigated methylation context-dependent links between genotype and gene expression in colorectal cancer (CRC) using paired tumor and normal-adjacent tissue (NAT) from 80 patients, thereby controlling for germline genomic background. By int…

</details>

---

### [Memory-safe high-performance sequence mapping with rammap](https://www.biorxiv.org/content/10.64898/2026.05.26.726289v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.32**
Published: 2026-05-29
Tags: sequence alignment, minimap2, software, Rust, pipelines

Presents a Rust reimplementation of minimap2 with performance and safety improvements and drop-in compatibility. Useful for infrastructure, but not specific to single-cell/omics interpretation, biomarkers, or translational oncology questions.

<details>
<summary>RSS summary</summary>

We introduce a reimplementation of the widely used mapping tool minimap2 in Rust called rammap. We demonstrate perfect concordance with minimap2, enabling its backwards compatibility as a drop-in replacement for minimap2-based workflows. Additionally, rammap implements performance optimizations for modern architectures and applications, including AVX512 and WASM v128 SIMD support for dynamic programming alignment and SIMD-accelerated chaining. These achieve comparable or better performance than …

</details>

---

### [Optical genome mapping identifies source-associated structural variant differences across early-passage human iPSCs](https://www.biorxiv.org/content/10.64898/2026.05.29.728843v1?rss=1)
*bioRxiv Genomics*  
Score: **0.30**
Published: 2026-05-31
Tags: structural variants, optical genome mapping, iPSC, genomics, model QC

Uses optical genome mapping to assess structural variants across iPSC sources, relevant to model quality control. It is not directly tied to oncology, therapy response, or single-cell computational methods in the title/summary provided.

<details>
<summary>RSS summary</summary>

Background: Induced pluripotent stem cells (iPSCs) are an important model for studying human diseases in vitro. However, previous studies have shown that iPSC reprogramming and extended cell culture can introduce genomic structural variants (SVs). Technologies like karyotyping, CNV microarrays, and whole-genome sequencing have limitations in resolution, sensitivity, or the ability to detect large and complex structural variants compared to optical genome mapping (OGM). OGM is a genome-wide struc…

</details>

---

### [QSyncFold: quantum neural network for multidimensional sync-discovery in protein folding](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag234/8698685?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.06**
Published: 2026-05-30
Tags: quantum ML, protein folding, method development

Proposes a hybrid quantum–classical neural network for protein structure prediction, which is far from single-cell/omics biomarker discovery and translational pediatric oncology workflows. The connection to the user’s stated priorities is not apparent from the title/abstract.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Quantum computing provides alternative encoding and sampling paradigms for protein structure prediction (PSP), but existing quantum-PSP methods are often limited by resource-scaling issues and by discrete or inefficient encodings for continuous coordinates. To address these limitations, we propose QSyncFold, a hybrid quantum–classical neural network framework that combines quantum superposition with differentiable learning. QSync…

</details>

---

## Other (9 shown / 9 total)

### [Histology-informed spatial domain identification through multi-view graph convolutional networks](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014281)
*PLOS Comp Bio*  
Score: **0.61**
Published: 2026-06-01
Tags: spatial transcriptomics, histology integration, graph convolution, clustering, tumor microenvironment

STESH integrates expression, spatial coordinates, and histology via multi-view graph convolution to identify spatial domains, aligning with spatial-omics computation relevant for tumor microenvironment studies. Not cancer-specific in the summary but the method is broadly applicable.

<details>
<summary>RSS summary</summary>

<p>by Huihui Zhang, Jiaxing Chang, Zirong Li, Yue Sun, Pinli Hu, Haoxiu Wang, Hang Yang, Yonglin Ren, Xingtan Zhang, Zehua Chen, Kok Wai Wong, Haojing Shao</p> Identifying spatial domains is crucial in spatial transcriptomics, yet effectively integrating gene expression, spatial location, and histology remains challenging. We present STESH, a Spatial Transcriptomics clustering method that combines Expression, Spatial information and Histology. STESH extracts histological features using a convolu…

</details>

---

### [Decoding Hierarchical Cell-Cell Communication in Spatial Multi-Omics with CellSTIC](https://www.biorxiv.org/content/10.64898/2026.05.27.728114v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.60**
Published: 2026-05-31
Tags: cell-cell communication, spatial multi-omics, ligand-receptor, tumor microenvironment, spatial transcriptomics

Focuses on inferring cell–cell communication while preserving spatial context using spatial transcriptomics/multi-omics, a common need in studying immunosuppressive niches and stromal interactions. The summary indicates hierarchical/region-specific signaling programs, which are transferable to oncology spatial analyses.

<details>
<summary>RSS summary</summary>

Cell-cell communication helps to coordinate tissue development, homeostasis, and immune responses, but identifying signaling interactions within intact tissues remains difficult. Although single-cell transcriptomics has enabled systematic inference of ligand-receptor interactions, dissociation disrupts spatial context and limits the identification of bona fide local signaling and region-specific communication programs. Spatial transcriptomics and spatial multi-omics offer the opportunity to stud…

</details>

---

### [Cophenetic Spatial Topology Embedding reveals multiscale tissue architecture in spatial omics](https://www.biorxiv.org/content/10.64898/2026.05.26.727847v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.55**
Published: 2026-05-29
Tags: spatial omics, tissue architecture, embedding, multiscale, computational method

Introduces a spatial omics embedding framework aimed at capturing multiscale tissue architecture without choosing a radius/neighborhood cutoff, which is a practical analysis advantage. Could be useful for summarizing tumor organization but no cancer-specific application is stated in the snippet.

<details>
<summary>RSS summary</summary>

The spatial organization of tissues emerges from cell interactions across multiple scales, yet current spatial omics analysis tools often emphasize local neighborhoods and may not summarize broader tissue architecture. Here we introduce Cophenetic Spatial Topology Embedding (COSTE), a computational framework that embeds directed nearest-neighbor distance profiles into a hierarchical metric space without requiring the user to define a spatial radius or neighborhood cutoff. COSTE can be applied to…

</details>

---

### [Supervised deep learning with gene functional annotation for cell classification](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014327)
*PLOS Comp Bio*  
Score: **0.52**
Published: 2026-06-01
Tags: scRNA-seq, cell classification, deep learning, interpretability, gene sets

Proposes supervised deep learning for scRNA-seq cell classification that incorporates gene functional annotation to improve interpretability beyond long DE gene lists. Useful for single-cell workflows, though not explicitly oncology-focused in the summary.

<details>
<summary>RSS summary</summary>

<p>by Zhexiao Lin, Yuanyuan Gao, Wei Sun</p> Gene-by-gene differential expression analysis is a widely used supervised approach for interpreting single-cell RNA-sequencing (scRNA-seq) data. However, modern scRNA-seq datasets often contain large numbers of cells, leading to the identification of many differentially expressed genes with extremely small p-values but negligible effect sizes, thus making biological interpretation difficult. To overcome this challenge, we developed Supervised Deep lea…

</details>

---

### [Epigenome Alterations and 3D Chromatin Architecture Remodeling in Inflammatory Macrophage Activation under Diabetic Conditions](https://www.biorxiv.org/content/10.64898/2026.05.27.727961v1?rss=1)
*bioRxiv Genomics*  
Score: **0.41**
Published: 2026-05-31
Tags: multi-omics, 3D genome, chromatin architecture, epigenomics, macrophages

Applies integrated multi-omics including 3D chromatin architecture to study macrophage inflammatory activation, showcasing analysis patterns relevant to chromatin/epigenome integration. However, it is not cancer-focused and is centered on diabetic conditions.

<details>
<summary>RSS summary</summary>

Aberrant monocyte and macrophage activation in diabetes drives chronic inflammation and complications. While epigenetic mechanisms are implicated, the role of 3D chromatin reorganization remains unclear. Using integrated multi-Omics, we profiled gene expression and 3D chromatin architecture in human CD14+ monocyte differentiated macrophages treated with high glucose plus TNF-alpha (HT) mimicking the diabetic milieu. HT induced inflammatory programs resembling those in diabetes, dynamically alter…

</details>

---

### [Genomics-Informed Approach Identifies Which Cell Types Regulate the Metabolome](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag330/8698028?rss=1)
*Bioinformatics (Oxford Academic)*  
Score: **0.40**
Published: 2026-05-29
Tags: metabolomics, mQTL, scRNA-seq atlas, data integration, cell-type mapping

Integrates large metabolite QTL resources with a large scRNA-seq atlas to map cell types regulating systemic metabolite levels, which is methodologically interesting for cross-dataset integration. Not oncology-focused and oriented toward body-wide metabolism rather than tumor context in the abstract snippet.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Metabolism occurs in a cell type-specific manner, but which cells regulate metabolite levels remains unclear. Here, we integrate some of the largest metabolite quantitative trait loci datasets, TOPMed and UK Biobank, with one of the most extensive single-cell RNA sequencing resources, Tabula Sapiens. This integration allows us to identify cell types that regulate metabolites body-wide. We find hepatocytes are the primary regulato…

</details>

---

### [Fibroblast-derived thrombospondin-1 shapes macrophage polarization in advanced human co-culture models](https://www.biorxiv.org/content/10.64898/2026.05.28.728363v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.32**
Published: 2026-05-29
Tags: tumor microenvironment, macrophage polarization, fibroblasts, organoid co-culture, mechanism

Describes advanced human co-culture/organoid-like models to study stromal–immune interactions and macrophage polarization, which is biologically relevant to TME. The snippet emphasizes model development/mechanism rather than computational single-cell/multi-omics analysis.

<details>
<summary>RSS summary</summary>

BackgroundTumor-associated macrophages (TAMs) are key drivers of the immunosuppressive tumor microenvironment (TME), supporting tumor progression through diverse functions. However, mechanistic studies of TAM polarization remain limited by the lack of physiologically relevant human model systems that capture stromal-immune interactions and macrophage heterogeneity. MethodsWe established advanced human co-culture systems that integrate healthy donor-derived macrophages with patient-derived organo…

</details>

---

### [Matched pancreatic cancer liver metastatic model system reveals cancer cell-dependent organotropism and site-specific tumor microenvironment reflective of human disease](https://www.biorxiv.org/content/10.64898/2026.05.27.728281v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.28**
Published: 2026-05-31
Tags: PDAC, metastasis model, organotropism, tumor microenvironment, preclinical model

Presents a matched primary/metastatic PDAC transplant model to study organotropism and site-specific TME, which is translational model-system work. Computational/omics components are not evident in the provided snippet.

<details>
<summary>RSS summary</summary>

Pancreatic ductal adenocarcinoma (PDAC) is a deadly, highly metastatic disease, driven by an interplay between cancer cells and the metastatic site-specific microenvironment. However, pre-clinical models that robustly capture these interactions within the context of matched primary and metastatic tumors are limited. Here, we present a novel transplant model system for matched pancreas and liver tumors to study PDAC metastatic progression. Using this model, we identified murine PDAC cell lines wi…

</details>

---

### [PCSK9 Exhibits Novel Nuclear Localization in LSEC and Its Targeting with Bioinspired Nanoparticles Reduces Colorectal Liver Metastasis](https://www.biorxiv.org/content/10.64898/2026.05.26.727886v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.20**
Published: 2026-05-29
Tags: metastasis, PCSK9, nanoparticles, liver microenvironment, colorectal cancer

Describes PCSK9 localization in liver sinusoidal endothelial cells and nanoparticle targeting to reduce colorectal liver metastasis. This is largely a mechanistic/nanoparticle delivery study with no clear computational or multi-omics component in the summary snippet.

<details>
<summary>RSS summary</summary>

Background & AimsColorectal cancer liver metastasis is the leading cause of mortality in affected patients, with liver sinusoidal endothelial cells playing a pivotal role in metastatic niche formation. Proprotein convertase subtilisin/kexin type 9 has emerged as a regulator of tumor biology, but its function in the hepatic microenvironment remains poorly defined. This study aimed to characterize the role and subcellular localization of PCSK9 in liver sinusoidal endothelial cells and to evaluate …

</details>

---
