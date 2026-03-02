# Weekly ToC Digest (week of 2026-03-02)

Triage focused on pediatric translational oncology, high-risk neuroblastoma relevance, and robust multi-omics/scRNA-seq methods. Direct NB papers are not in this set; top items emphasize translational cancer biology, therapy resistance mechanisms, and transferable computational frameworks for multi-omics and single-cell data. Prioritized pediatric NB-related and translational multi-omics items; highlighted NB biomarker study and several translational computational-method papers with potential clinical utility. Triage of weekly TOC items for pediatric translational oncology and computational biology. NB-focused items are limited in this set; prioritized entries emphasize single-cell and multi-omics methods with translational or interpretability value, plus biomarker and therapy-resistance insights that could transfer to pediatric contexts. NB-specific papers are scarce this week. Highlighted items are computational/chromatin-regulatory studies with potential translational relevance (imm

**Included:** 20 (score ≥ 0.35)  
**Scored:** 30 total items

**Models:** `gpt-5-nano` first-pass, `gpt-5.2` re-rank top 40

---

## Neuroblastoma (2 shown / 2 total)

### [Spatially Integrated Multi-Omics reveals the Multicellular Landscape of Progenitor-Driven Glioblastoma Progression](https://www.biorxiv.org/content/10.64898/2026.02.26.708154v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.71**  
Published: 2026-02-28T00:00:00+00:00
Tags: multi-omics, spatial transcriptomics, single-cell, tumor microenvironment, progression, computational biology

Uses spatially integrated multi-omics to connect genomic drivers, cell states, and microenvironment interactions with progression/survival in glioblastoma. While not pediatric/NB-specific, the integration of spatial + omics for multicellular programs is highly transferable to tumor ecosystem and biomarker analyses.

<details>
<summary>RSS summary</summary>

Glioblastoma is the most lethal primary brain tumor, driven by complex interactions between plastic malignant cells and a diverse tumor microenvironment. Despite advances in single-cell profiling, how genomic drivers and the tumor microenvironment interact to determine tumor progression and patient survival remains poorly understood. While cellular states have been cataloged, the multicellular logic coordinating these into lethal phenotypes remains unresolved. Here, we integrate whole-exome sequ…

</details>

---

### [Pembrolizumab in advanced malignant peripheral nerve sheath tumors: a single-arm phase 2 trial](https://www.nature.com/articles/s41698-026-01303-6)
*npj Precision Oncology*  
Score: **0.48**  
Published: 2026-02-27T00:00:00+00:00
Tags: clinical trial, pembrolizumab, immunotherapy, translational, MPNST

A phase 2 immunotherapy clinical trial is translationally relevant and may include response/biomarker correlates, but the TOC text provides no molecular/computational details and the disease context is not neuroblastoma. Scored as moderate due to clinical response relevance only.

<details>
<summary>RSS summary</summary>

<p>npj Precision Oncology, Published online: 27 February 2026; <a href="https://www.nature.com/articles/s41698-026-01303-6">doi:10.1038/s41698-026-01303-6</a></p>Pembrolizumab in advanced malignant peripheral nerve sheath tumors: a single-arm phase 2 trial

</details>

---

## AI (5 shown / 7 total)

### [MAP: A Knowledge-driven Framework for Predicting Single-cell Responses for Unprofiled Drugs](https://www.biorxiv.org/content/10.64898/2026.02.25.708091v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.86**  
Published: 2026-02-27T00:00:00+00:00
Tags: single-cell, perturbation modeling, drug response, zero-shot, knowledge graph, AI/ML, translational

Proposes a knowledge-integrated perturbation model aimed at zero-shot prediction of single-cell responses for unprofiled drugs, directly aligned with drug response/resistance modeling workflows. The summary emphasizes encoding mechanistic relationships rather than treating drugs as IDs, supporting transfer to translational perturbation settings.

<details>
<summary>RSS summary</summary>

Predicting how cells respond to chemical perturbations is one of the goals for building virtual cells, yet experimentally profiled compounds cover only a small fraction of this space. Existing models struggle to generalize to unprofiled compounds, as they typically treat drugs as isolated identifiers without encoding their mechanistic relationships. We present MAP, a framework that integrates structured biological knowledge into cellular perturbation modeling and supports zero-shot prediction fo…

</details>

---

### [Uncertainty-aware synthetic lethality prediction with pretrained foundation models](https://www.biorxiv.org/content/10.64898/2026.02.25.708096v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.78**  
Published: 2026-02-27T00:00:00+00:00
Tags: synthetic lethality, foundation model, uncertainty, AI/ML, target prioritization, combination therapy

Introduces a synthetic-lethality prediction framework leveraging pretrained biological foundation models with uncertainty awareness, relevant to prioritizing combination targets in cancer. The focus on generalization beyond curated PPI/GO can be useful for discovery in less-studied pediatric contexts.

<details>
<summary>RSS summary</summary>

Synthetic lethality (SL) offers a promising paradigm for targeted cancer therapy, yet experimental identification of SL gene pairs remains costly, context-dependent, and biased toward well-studied genes. Existing computational approaches often rely on curated protein-protein interaction (PPI) networks and Gene Ontology (GO) annotations, which limit their ability to generalize to novel genes. Here we introduce CILANTRO-SL, a two-stage, graph-free framework that leverages pretrained biological fou…

</details>

---

### [mSWI/SNF complex inhibition sensitizes KRAS-mutant lung cancers to targeted therapies via epithelial-mesenchymal subversion](https://www.biorxiv.org/content/10.64898/2026.02.27.708377v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.54**  
Published: 2026-03-01T00:00:00+00:00
Tags: therapy resistance, EMT, cell-state plasticity, chromatin remodeling, combination therapy, targeted therapy

Mechanistic/translational resistance study identifying chromatin remodeling (mSWI/SNF) as a determinant of EMT-mediated targeted-therapy resistance and testing a clinical-grade inhibitor to dampen resistance. Cancer type differs, but EMT/cell-state plasticity framing is relevant to resistance programs broadly.

<details>
<summary>RSS summary</summary>

Targeted therapies for KRAS-mutant non-small lung cancer (NSCLC) have shown promising clinical results, however, incomplete tumoral responses and the inevitable emergence of therapeutic resistance remain critical challenges. Here we identify mSWI/SNF chromatin remodeling complexes as critical determinants of (EMT)-mediated KRAS inhibitor inefficacy and resistance in KRAS G12C lung cancers. Treatment with the clinical-grade SMARCA4/2 inhibitor, FHD-286, dampens EMT-mediated acquired resistance in…

</details>

---

### [Tumor-Infiltrating lymphocyte dynamics as biomarkers of neoadjuvant treatment response in luminal breast cancer](https://link.springer.com/article/10.1186/s12885-026-15795-9)
*BMC Cancer*  
Score: **0.40**  
Published: 2026-02-27T00:00:00+00:00
Tags: biomarkers, treatment response, TILs, immunology, breast cancer

Focuses on TIL dynamics as biomarkers of neoadjuvant response, which is conceptually aligned with response prediction. The RSS entry lacks a summary, so the computational rigor and biomarker validation details are unclear.

---

### [Expanding the RB1 variant landscape of heritable retinoblastoma: unlocking precision oncology potential in Southern Africa](https://link.springer.com/article/10.1186/s12885-026-15789-7)
*BMC Cancer*  
Score: **0.38**  
Published: 2026-02-28T00:00:00+00:00
Tags: pediatric oncology, retinoblastoma, RB1, germline variants, precision oncology

Pediatric oncology adjacent (retinoblastoma) and focused on germline RB1 variant landscape, relevant to precision oncology and variant interpretation. No summary is provided here, so it is unclear how much computational analysis or clinical actionability is included.

---

## Methods (5 shown / 11 total)

### [A unified framework for multiomics deconvolution](https://www.nature.com/articles/s41592-026-03008-x)
*Nature Methods*  
Score: **0.92**  
Published: 2026-03-02T00:00:00+00:00
Tags: multi-omics, deconvolution, translational, Nature Methods

Introduces a universal multi-omics deconvolution framework; directly relevant to biomarker discovery and mechanism inference across omics.

---

### [DECODE: deep learning-based common deconvolution framework for various omics data](https://www.nature.com/articles/s41592-026-03007-y)
*Nature Methods*  
Score: **0.90**  
Published: 2026-03-02T00:00:00+00:00
Tags: deconvolution, multi-omics, deep learning, translational

DL-based universal deconvolution across transcriptomics, proteomics, and metabolomics; aligns with multi-omics biomarker and pathway activity analysis.

---

### [scDock: Streamlining drug discovery targeting cell–cell communication via scRNA-seq analysis and molecular docking](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag103/8503425?rss=1)
*Bioinformatics (Oxford Academic)*  
Score: **0.82**  
Published: 2026-03-02T00:00:00+00:00
Tags: scRNA-seq, cell-cell communication, drug discovery, molecular docking, pipeline, translational

Presents an integrated pipeline connecting scRNA-seq-based cell–cell communication analysis to structure-based molecular docking, directly bridging single-cell computation to drug discovery. This is methodologically transferable for identifying actionable ligand–receptor targets and prioritizing compounds.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle">Summary</div>Identifying drugs that target intercellular communication networks represents a promising therapeutic strategy, yet linking single-cell RNA sequencing (scRNA-seq) analysis to structure-based drug screening remains technically challenging and requires substantial bioinformatics expertise. We present scDock, an integrated and user-friendly pipeline that seamlessly connects scRNA-seq data processin…

</details>

---

### [Achieving spatial multi-omics integration from unaligned serial sections with DIME](https://www.biorxiv.org/content/10.64898/2026.02.25.707961v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.78**  
Published: 2026-02-28T00:00:00+00:00
Tags: spatial omics, multi-omics integration, representation learning, embeddings, computational method

Introduces DIME for "diagonal integration" of spatial multi-omics collected from unaligned serial sections, addressing a common real-world setting where modalities do not share features. The focus on integrated embeddings is broadly useful for multi-omics representation learning and downstream biomarker discovery workflows.

<details>
<summary>RSS summary</summary>

Learning integrated representations from spatial multi-omics data is a fundamental challenge, particularly in the context of "diagonal integration", where data are collected from serial tissue sections across distinct omics modalities. Existing methods typically rely on the assumption of feature intersection to construct a common metric space, a prerequisite that is absent in this setting. To address this, we propose the Diagonal Integration Model for Spatial Multi-omics Embedding (DIME), a nove…

</details>

---

### [Hypoxia and Associated Acidosis Generate Cell-Type Specific Myeloid Responses in Glioblastoma](https://www.biorxiv.org/content/10.64898/2026.02.26.707379v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.67**  
Published: 2026-02-28T00:00:00+00:00
Tags: scRNA-seq, spatial transcriptomics, tumor microenvironment, myeloid cells, hypoxia, multi-omics

Integrates cyclic IHC, scRNA-seq, spatial transcriptomics, in vitro models, and DNA methylation profiling across many GBMs to dissect hypoxia/acidosis effects on myeloid states. The multi-modal design and microenvironmental-state readouts are transferable to pediatric solid tumor TME and response studies.

<details>
<summary>RSS summary</summary>

Hypoxia is a defining feature of glioblastoma (GBM), yet how it cooperates with hypoxia-associated acidosis to shape microglia and infiltrating monocyte-derived macrophages (MDM) remains poorly understood. We integrated cyclic immunohistochemistry, single-cell RNA sequencing, spatial transcriptomics, in vitro cell cultures, and DNA methylation profiling to outline hypoxia-driven responses in up to 136 GBMs. These hypoxic niches were selectively enriched for MDMs that activated carbonic anhydrase…

</details>

---

## Other (0 shown / 0 total)
