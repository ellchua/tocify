# Weekly ToC Digest (week of 2026-11-02)

This digest is automatically generated from this week's RSS items and categorized into Neuroblastoma, AI, Methods, and Other.

**Included:** 26 (all ranked papers)  
**Scored:** 26 total items

**Models:** `gpt-5-nano` first-pass, `gpt-5.2` re-rank top 40

---

## Neuroblastoma (14 shown / 14 total)

### [Cross-modal transfer enables zero-shot prediction of drug responses across cell types](https://www.biorxiv.org/content/10.64898/2026.10.27.762107v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.84**
Published: 2026-11-01
Tags: drug response, perturbation, transfer learning, cross-modal, single-cell, computational, therapy response

Introduces GRAFT to transfer information between genetic and chemical perturbation modalities for zero-shot drug-response prediction across cellular contexts. Directly relevant to modeling therapy response/resistance using perturbation data and transferable regulatory signals.

<details>
<summary>RSS summary</summary>

Genetic perturbation responses capture the functional consequences of gene-level interventions, providing functional information for drug-response prediction beyond baseline gene expression. Within the same cellular context, genetic and chemical perturbations engage shared regulatory networks, providing a biological basis for transferring response information between perturbation modalities. Here we introduce GRAFT, a framework that combines genetic-to-chemical transfer with cross-context predic…

</details>

---

### [Single-cell perturbation prediction near experimental reproducibility with GENESIS](https://www.biorxiv.org/content/10.64898/2026.10.30.763081v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.82**
Published: 2026-11-01
Tags: single-cell, perturbation, generative model, dose response, drug response, computational

Presents GENESIS, a generative energy-based framework predicting single-cell response distributions from experimental covariates (context, perturbation, dose). Highly relevant to scalable perturbation modeling for mechanism and drug-response prediction.

<details>
<summary>RSS summary</summary>

Single-cell perturbation atlases now measure cellular responses across millions of cells, but the number of possible combinations of biological contexts, perturbations, doses and other experimental variables make exhaustive measurements infeasible. Here, we present GENESIS (Generative ENErgy-based Single-cell Intervention Simulator), a generative framework that predicts single-cell response distributions directly from the covariates defining an experimental condition. GENESIS uses a lightweight …

</details>

---

### [Persistent Malignant Stem Cells Drive Relapse After Transplantation in the Myelodysplastic Neoplasms](https://www.biorxiv.org/content/10.64898/2026.10.28.761431v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.74**
Published: 2026-10-29
Tags: biomarkers, relapse prediction, minimal residual disease, ultra-low-input sequencing, clinical cohort, translational

Develops an ultra-low-input DNA sequencing assay applied to a multicenter post-transplant cohort, where mutation detection in HSCs preceded relapse and strongly predicted relapse risk. This is a clinically oriented minimal residual disease/relapse-biomarker framework that is transferable to pediatric oncology settings even though it is not neuroblastoma-focused.

<details>
<summary>RSS summary</summary>

The myelodysplastic neoplasms (MDS) are initiated in hematopoietic stem cells (HSCs). Allogeneic transplantation is the only curative treatment, but relapse remains common. To detect rare MDS HSCs, we developed an ultra-low-input DNA sequencing assay and applied it to a multicenter post-transplant cohort. HSC mutation detection preceded relapse by a mean of 98 days and strongly predicted relapse (95.7% vs. 30.3%; hazard ratio, 8.8; P=1.5e-7). Mutations detected only in unfractionated marrow did …

</details>

---

### [SpatioCAD: Context-aware graph diffusion model for pinpointing spatially variable genes in heterogeneous tissues](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014923)
*PLOS Comp Bio*  
Score: **0.66**
Published: 2026-11-02
Tags: spatial transcriptomics, spatially variable genes, graph diffusion, tumor heterogeneity, biomarkers, computational

Presents a context-aware graph diffusion model for detecting spatially variable genes while addressing confounding from cell density variations in heterogeneous tissues such as tumors. Transferable to tumor spatial marker discovery and region-specific programs.

<details>
<summary>RSS summary</summary>

<p>by Hao Wen, Qunlun Shen, Shuqin Zhang</p> Spatial transcriptomics enables comprehensive characterization of tissue architecture, and the identification of spatially variable genes (SVGs) is a critical step for defining region-specific molecular markers and uncovering spatially regulated mechanisms across diverse biological contexts. However, most existing methods for SVG detection overlook cell density variations, a major confounding factor in complex tissues such as tumors, where heterogeneo…

</details>

---

### [Somatic variant detection using an individual-specific pangenome reference: a case study of a pediatric acute myeloid leukemia patient](https://www.biorxiv.org/content/10.64898/2026.10.27.756267v1?rss=1)
*bioRxiv Genomics*  
Score: **0.60**
Published: 2026-10-30
Tags: pediatric AML, somatic variants, pangenome, Nanopore, single-cell, clinical genomics

Evaluates somatic variant calling using an individual-specific pangenome built from Nanopore assemblies in a pediatric AML case, including bulk and single-cell contexts. Relevant to pediatric cancer genomics pipelines and reducing reference bias in clinically relevant genes.

<details>
<summary>RSS summary</summary>

Background: Reference bias can cause false or missed somatic variant calls, also in clinically relevant genes. Pangenomes aim to reduce this bias, but evidence of their diagnostic benefit is lacking. We tested how an individual-specific pangenome affects somatic variant calling in bulk samples and single cells. Results: We built an individual-specific pangenome from Nanopore assemblies of normal hematopoietic stem and progenitor cells of a child with acute myeloid leukemia. We aligned short-read…

</details>

---

### [Mapping chromothripsis across human cancer cell lines reveals genome evolution and therapeutic vulnerabilities](https://www.biorxiv.org/content/10.64898/2026.10.30.763015v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.54**
Published: 2026-10-30
Tags: chromothripsis, copy number, cancer cell lines, pipeline, therapeutic vulnerabilities

Develops a pipeline to detect chromothripsis-like events from copy-number data and maps them across 1,646 cancer cell lines, linking patterns to therapeutic vulnerabilities. Method and vulnerability mapping are broadly relevant, but not pediatric/NB-specific per the summary.

<details>
<summary>RSS summary</summary>

Chromothripsis is a catastrophic chromosomal rearrangement found in ~40% of primary tumors, yet its landscape in cancer cell lines (CCLs) remains uncharacterized. To address this gap, we developed a pipeline to detect chromothripsis-like events from copy-number data alone, and mapped their occurrence across 1,646 human CCLs. Chromothripsis patterns were broadly conserved between primary tumors and CCLs, but CCLs displayed some distinct features and a significantly higher prevalence of chromothri…

</details>

---

### [IsoPrism predicts isoform usage from short-read transcriptomes through structure-guided competitive allocation](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag621/9006417?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.52**
Published: 2026-11-02
Tags: isoforms, short-read, single-cell, spatial transcriptomics, transcript structure, computational

Introduces IsoPrism to infer isoform proportions from gene-level short-read single-cell/spatial expression using transcript structural annotations. Potentially useful for isoform-level biomarker hypotheses when long-read is unavailable, though not cancer-specific.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Short-read (SR) single-cell and spatial transcriptomics characterize gene expression but provide limited resolution of within-gene isoform usage. Here we present IsoPrism, a structure-guided competitive allocation framework that predicts isoform proportions from gene-level expression. IsoPrism combines individual expression and local expression context with transcript structural annotations, using shared iterative updates to pres…

</details>

---

### [ERX-208 Induces Endoplasmic Reticulum Stress and Suppresses DNA Damage Repair to Enhance Chemotherapy Sensitivity](https://aacrjournals.org/cancerrescommun/article/6/10/2541/789712/ERX-208-Induces-Endoplasmic-Reticulum-Stress-and)
*AACR Cancer Research Communications*  
Score: **0.44**
Published: 2026-10-29
Tags: DNA damage repair, chemotherapy sensitization, ER stress, combination therapy, resistance

Mechanistic therapy-sensitization study linking ER stress induction and suppression of DNA damage repair to improved chemotherapy sensitivity, but in ovarian cancer models. Relevant mainly as an adjacent DDR/combination-therapy concept rather than pediatric/neuroblastoma-specific translational work.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle"></div>Ovarian cancer remains the most lethal gynecologic malignancy, driven in part by intrinsic and acquired therapeutic resistance to DNA-damaging chemotherapies. We previously demonstrated that ERX-208 targets the protein encoded by lysosomal acid lipase to induce catastrophic endoplasmic reticulum (ER) stress in ovarian cancer models. As ERX-class compounds advance to clinical development, rational combi…

</details>

---

### [PRIDA: Cross-domain protein-RNA interaction prediction based on biological language models and conditional domain-adversarial learning](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag802/9004413?rss=1)
*Bioinformatics (Oxford Academic)*  
Score: **0.33**
Published: 2026-10-29
Tags: biological language models, protein-RNA interactions, domain adaptation, machine learning

Introduces a biological language model approach for protein–RNA interaction prediction with domain-adversarial learning, targeting generalizability across domains. It is computationally interesting but not clearly tied to single-cell/omics biomarker discovery or translational oncology use cases from the title/abstract snippet.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle">Motivation</div>The interaction between RNA and proteins is a complex process that involves atom-level details and dynamic changes. Accurate prediction of protein-RNA interactions helps researchers to develop protein-targeted small nucleic acid drugs. However, the scarcity of protein-RNA complexes and experimental data limits the ability of existing data-driven supervised learning methods to learn latent pat…

</details>

---

### [Myogenic and neural transcriptional circuitry determines PAX3::FOXO1 driven initiation of rhabdomyosarcoma](https://www.biorxiv.org/content/10.64898/2026.10.27.761987v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.30**
Published: 2026-10-29
Tags: rhabdomyosarcoma, PAX3-FOXO1, pediatric cancer, cell identity, model system

Uses an inducible/degradable PAX3::FOXO1 iPSC-derived ARMS model to study initiation and cell identity. Pediatric solid tumor relevance, but the summary emphasizes a model/mechanism study and does not mention omics integration or biomarker/response modeling.

<details>
<summary>RSS summary</summary>

Alveolar rhabdomyosarcoma (ARMS) is largely characterized by the PAX3::FOXO1 (P3F) fusion oncoprotein, yet the mechanisms by which P3F initiates tumorigenesis and cell identity remain incompletely understood. Here, we describe a doxycycline inducible, degradable PAX3::FOXO1 iPSC derived ARMS model (ddP3F) to interrogate the temporal requirements of P3F for ARMS cell identity and proliferation. We show that, while P3F is essential for initiating ARMS transformation from precursor cells, continued…

</details>

---

### [AxialTAD: an axial contrast architecture for TAD boundary prediction from Hi-C](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag603/9006411?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.30**
Published: 2026-11-02
Tags: Hi-C, TAD boundaries, deep learning, 3D genome

Proposes a compact deep learning model for predicting TAD boundaries from Hi-C matrices. Methodologically interesting but not clearly tied to translational oncology or omics biomarker workflows in the provided abstract snippet.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Topologically associating domains (TADs) are fundamental units of 3D chromatin organization, and accurate identification of TAD boundaries from Hi-C contact matrices remains an active challenge. Recent deep learning approaches have advanced substantially, but these rely on increasingly large architectures, raising the question of how much of this capacity boundary-centric analysis requires. Here we propose AxialTAD, a compact dee…

</details>

---

### [APR-246 and DNA polymerase theta inhibitor Novobiocinsynergize to eliminate TP53 mutant AML cells](https://www.biorxiv.org/content/10.64898/2026.10.30.763043v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.28**
Published: 2026-11-02
Tags: AML, TP53, combination therapy, APR-246, POLQ inhibitor

Reports a combination strategy (APR-246 with POLθ inhibitor novobiocin) against TP53-mutant AML cells. Translationally themed, but leukemia- and wet-lab–leaning from the summary with no clear computational/omics component described.

<details>
<summary>RSS summary</summary>

Patients with TP53-mutant myelodysplasic syndromes (MDS) and acute myeloid leukemia (AML) are usually elderly and frail with highly resistant malignant clones and a median survival of less than one year. The molecule APR-246 targets mutant p53 but also perturbs redox homeostasis. APR-246 showed favorable safety in a phase 3 trial in TP53 mutant MDS/AML but did not meet efficacy endpoints in combination with azacitidine (AZA), highlighting the need for more effective combination strategies. Conve…

</details>

---

### [Identification of TM4SF1 as a novel tumor-associated antigen in biliary tract cancers targetable by CAR-T cell therapy](https://www.biorxiv.org/content/10.64898/2026.10.28.762675v1?rss=1)
*bioRxiv Cancer Biology*  
Score: **0.18**
Published: 2026-10-30
Tags: tumor antigen, CAR-T, biliary tract cancer, immunotherapy

Identifies TM4SF1 as a tumor-associated antigen in biliary tract cancers with CAR-T targeting rationale. Solid tumor immunotherapy angle but outside pediatric/NB scope and no computational/omics emphasis is evident from the summary snippet.

<details>
<summary>RSS summary</summary>

Background & Aims: Biliary tract cancers (BTCs) are a rare set of genetically heterogenous and aggressive malignancies associated with delayed presentation and poor prognosis. Thus, there is an unmet clinical need for the development of novel therapeutic strategies. Identification of uniquely expressed cell-surface, tumor-associated antigens (TAAs) holds promise in epithelial tumors more generally, as they can serve as ligands for a variety of therapeutics, including CAR-T cell therapy. However,…

</details>

---

### [Cd47 deficiency in mouse brain under excitotoxic stress](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014897)
*PLOS Comp Bio*  
Score: **0.08**
Published: 2026-10-29
Tags: snRNA-seq, mouse brain, CD47, neuroscience

Single-nucleus RNA-seq study of Cd47 deficiency in mouse brain under excitotoxic stress. Not oncology-focused and limited translational relevance to pediatric cancer from the title/summary.

<details>
<summary>RSS summary</summary>

<p>by Mehrdad Kiamari, Amalia Perna, Matteo Santoro, Abdel El Ariss, Eloise Berson, Jing Zhao, Chandresh R. Gajera, Edward J. Fox, Kathleen S. Montine, Nima Aghaeepour, Thomas J. Montine</p> CD47, a transmembrane protein with well-established immune and neurodevelopmental roles, was recently identified as an abundant presynaptic protein in humans, non-human primates, and mice. To investigate how <i>Cd47</i> influences neuronal stress and injury, we performed single-nucleus RNA sequencing of hipp…

</details>

---

## AI (0 shown / 0 total)

## Methods (2 shown / 2 total)

### [Verix: a toolkit for benchmarking and harmonization of complex structural variants](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag803/9004412?rss=1)
*Bioinformatics (Oxford Academic)*  
Score: **0.63**
Published: 2026-10-29
Tags: structural variants, benchmarking, toolkit, genomics, reproducible pipelines

Presents a benchmarking/harmonization toolkit targeted to complex structural variants, addressing evaluation gaps beyond standard SV frameworks. This is methodologically relevant for genomic biomarker discovery pipelines where SVs/complex rearrangements may impact oncogenic drivers or resistance mechanisms.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle">Summary</div>Complex structural variants (CSVs) represent a critical yet under-characterized class of large-scale genomic rearrangements. CSV discovery requires not only the accurate detection of individual breakpoints but also their assembly into a single, structurally resolved CSV event. However, CSV method development is currently limited by a diagnostic gap: existing SV evaluation frameworks are designed…

</details>

---

### [HiCP2GAN: a plug and play foundation model-based GAN for Hi-C enhancement](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag610/9006419?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.34**
Published: 2026-11-02
Tags: Hi-C, GAN, 3D genomics, deep learning, data enhancement

Describes a foundation model-based GAN approach to enhance low-resolution Hi-C contact maps. Primarily a 3D genomics methods paper with unclear immediate connection to single-cell tumor biology or biomarker discovery from the summary.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>The 3D organization of chromatin shapes gene regulation and cellular function, and Hi-C has become the primary technique for mapping chromatin interactions genome-wide. However, high-resolution data remain costly and scarce, leaving many studies with sparse contact maps that limit downstream analysis. Deep learning methods, especially generative adversarial networks (GANs), have shown promise for enhancing low-resolution Hi-C dat…

</details>

---

## Other (10 shown / 10 total)

### [Building a pre-trained model for single-cell multi-omics analysis beyond highly variable feature selection](https://www.nature.com/articles/s41467-026-78933-8)
*Nature Communications*  
Score: **0.90**
Published: 2026-11-02
Tags: single-cell, multi-omics, pretrained model, scFM, representation learning, integration, computational

Presents scPMMI, a pre-trained model for single-cell multi-omics that explicitly aims to capture signal beyond highly variable feature selection and supports multiple downstream tasks. Strong match to single-cell foundation/representation learning and multi-omics integration needs.

<details>
<summary>RSS summary</summary>

<p>Nature Communications, Published online: 02 November 2026; <a href="https://www.nature.com/articles/s41467-026-78933-8">doi:10.1038/s41467-026-78933-8</a></p>The authors present scPMMI, a pre-trained model that learns from multiple types of molecular data in individual cells. They show that this method captures signals beyond highly variable features and supports diverse multi-omics analysis tasks.

</details>

---

### [EMBER: Contextualized gene representations enable compositional prediction across cellular perturbations](https://www.biorxiv.org/content/10.64898/2026.10.28.762655v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.80**
Published: 2026-10-31
Tags: foundation model, gene representations, perturbation, combinatorial, representation learning, single-cell

Describes EMBER, a gene-centric foundation model learning context-dependent gene representations to predict effects of single and combined perturbations. Useful for compositional perturbation modeling and identifying context-specific vulnerabilities.

<details>
<summary>RSS summary</summary>

Predicting the effects of individual and combined perturbations is central to identifying context-specific vulnerabilities and designing targeted interventions. Yet learning gene representations that integrate functional relationships with cellular state and transfer across perturbation tasks remains challenging. Here, we present EMBER, a gene-centric foundation model that learns context-dependent functional representations of genes. EMBER combines interaction-guided contrastive learning with ma…

</details>

---

### [Multi-omics survival prediction under modality incompleteness via dynamic modality graphs and product-of-experts fusion](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag596/9006413?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.74**
Published: 2026-11-02
Tags: multi-omics, survival prediction, missing data, data integration, clinical modeling, machine learning

Proposes DMG-PoE for multi-omics survival prediction explicitly handling modality missingness via a dynamic modality graph and product-of-experts fusion. Strongly transferable to clinical/translational cohort modeling where incomplete omics is common.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Multi-omics data integration provides rich complementary information for cancer survival prediction, yet modality incompleteness is pervasive in real-world medical data and undermines the stability and applicability of existing models. To address this challenge, we propose DMG-PoE, a robust learning framework for multi-omics survival prediction. The framework introduces a dynamic modality graph (DMG) that propagates complementary…

</details>

---

### [SpaceBlocks: an expert-supervised workflow for cohort-level spatial transcriptomics analysis](https://www.biorxiv.org/content/10.64898/2026.10.30.763100v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.70**
Published: 2026-11-02
Tags: spatial transcriptomics, workflow, Snakemake, reproducibility, cohort-level, computational pipeline

Introduces SpaceBlocks, a Snakemake workflow for reproducible and scalable cohort-level spatial transcriptomics with expert-supervised decision points. High practical value for robust pipelines and provenance in translational tissue studies.

<details>
<summary>RSS summary</summary>

Spatial transcriptomics (ST) has become a central tool in tissue biology, with multiple technologies addressing the issue from a different angle. Each platform imposes its own raw data format conventions, and few workflows address reproducible, scalable, cohort-level analyses while integrating expert-supervised decisions to strengthen computational and biological provenance. In this work, we present SpaceBlocks, a Snakemake workflow for reproducible cohort-level analysis of single-cell resolutio…

</details>

---

### [Enhancing cell type deconvolution in spatial transcriptomics through incorporating multiple reference panels](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1012345)
*PLOS Genetics*  
Score: **0.64**
Published: 2026-11-02
Tags: spatial transcriptomics, deconvolution, multi-reference, single-cell reference, tumor microenvironment

Focuses on improving spatial transcriptomics cell-type deconvolution by incorporating multiple scRNA-seq reference panels to handle reference mismatch. Relevant to robust tumor microenvironment and cell-state quantification in clinical specimens.

<details>
<summary>RSS summary</summary>

<p>by Peiyao Zhao, Justin A. Colacino, Ying Ma, Xiang Zhou</p> Spatial transcriptomics (ST) captures gene expression with spatial context but often at resolutions where each spatial location contains multiple cells. Deconvolution of cell-type proportions is therefore essential, yet existing methods typically rely on a single scRNA-seq reference, limiting accuracy when references are incomplete or mismatched especially in complex tissues. These limitations highlight the need for a framework that …

</details>

---

### [Single-cell foundation models transfer across deep evolutionary distances](https://www.biorxiv.org/content/10.64898/2026.10.28.762457v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.55**
Published: 2026-10-31
Tags: single-cell, foundation models, zero-shot, embeddings, benchmarking

Assesses zero-shot embeddings from Geneformer, scGPT and UCE on distant species developmental atlases without retraining. Useful for understanding scFM generalization, though not directly tied to oncology or clinical endpoints from the summary.

<details>
<summary>RSS summary</summary>

Single-cell foundation models are trained mostly on human and other vertebrate cells, yet are proposed as general-purpose models. Whether they transfer far beyond their training species is unknown. We evaluated zero-shot embeddings from Geneformer, scGPT and UCE against PCA and scVI in developmental atlases of Dictyostelium discoideum and Caenorhabditis elegans, separated from humans by 1.2 billion and 600 million years. Without retraining, the models recovered developmental order and preserved …

</details>

---

### [Transcript quantification using context-aware MAPQ modeling](https://www.biorxiv.org/content/10.64898/2026.10.26.761703v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.43**
Published: 2026-10-31
Tags: long-read RNA-seq, transcript quantification, MAPQ, statistical modeling, pipeline

Presents CALMeR, a context-aware model for long-read transcript quantification that avoids static MAPQ filtering by conditioning on transcript/read properties. Useful sequencing-methods improvement, but indirect to oncology/translational questions.

<details>
<summary>RSS summary</summary>

Long-read RNA-sequencing enables accurate transcript-level quantification, but prior static filtering of MAPping Quality (MAPQ) can ignore transcriptome-intrinsic ambiguity, and discarding reads on this basis distorts abundance estimates. We present CALMeR (Context-Aware Long-read MAPQ-based Rescuer), which conditions MAPQ filtering on transcript properties rather than a fixed cutoff. CALMeR fits a zero-one inflated beta regression on various transcript and read-level features to model the expec…

</details>

---

### [ADInteractome: discovering Alzheimer’s disease-associated gene interactions from single-nucleus RNA-seq data through neural interaction states](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag601/9006415?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.22**
Published: 2026-11-02
Tags: snRNA-seq, gene interactions, interpretable AI, network inference, Alzheimer's

Introduces an interpretable AI framework to infer gene interactions from snRNA-seq in Alzheimer’s disease. Method may be transferable, but the application domain is neurodegeneration rather than cancer/translational oncology.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Alzheimer’s disease (AD) is a prevalent neurodegenerative disorder characterized by amyloid $\beta $ plaques and tau tangles, yet its molecular mechanisms remain elusive. Advances in single-nucleus RNA sequencing (snRNA-seq) have provided profound insights into the cellular heterogeneity of AD-associated brains. In this study, we introduce ADInteractome, an interpretable artificial intelligence framework designed to discover gene…

</details>

---

### [Ten simple rules for teaching practical AI to health data teams](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014915)
*PLOS Comp Bio*  
Score: **0.22**
Published: 2026-10-30
Tags: AI in healthcare, best practices, training, governance

General guidance piece on operationalizing AI in healthcare data teams (privacy, compliance, expectations) rather than a method for omics, single-cell analysis, or translational oncology discovery. Useful for team/process considerations but low direct project relevance.

<details>
<summary>RSS summary</summary>

<p>by Eric Just, Aik Choon Tan</p> Artificial intelligence (AI) has rapidly entered the healthcare discourse, often framed through the lens of frontier models, generative systems, and transformative potential. Yet for many healthcare data teams, the practical path from curiosity to capability remains unclear. Concerns about data privacy, regulatory compliance, brittle systems, and misaligned expectations frequently stall progress before meaningful experimentation can begin. This paper presents t…

</details>

---

### [From 16S rDNA-seq count data to microbial interactions: A comprehensive benchmark of co-occurrence network inference methods](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014898)
*PLOS Comp Bio*  
Score: **0.12**
Published: 2026-10-29
Tags: microbiome, co-occurrence networks, benchmarking, 16S sequencing

Benchmarks co-occurrence network inference methods for microbial 16S count data; the application domain is microbiome ecology rather than cancer or single-cell tumor biology. Network benchmarking principles may be broadly informative but are not closely aligned with the stated oncology priorities.

<details>
<summary>RSS summary</summary>

<p>by Piero Mariotto, Matteo Baldan, Barbara Di Camillo, Giacomo Baruzzo</p> The study of the ecological interactions within microbial communities is crucial for gaining deeper insights into community ensemble and evolution, enabling targeted modification of the community’s functionality. Many inference approaches have been developed to derive ecological associations. Here, we focus on co-occurrence computational methods, in which associations are inferred from metagenomic abundance data by meas…

</details>

---
