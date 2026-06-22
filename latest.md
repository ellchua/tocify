# Weekly ToC Digest (week of 2026-06-22)

This digest is automatically generated from this week's RSS items and categorized into Neuroblastoma, AI, Methods, and Other.

**Included:** 36 (all ranked papers)  
**Scored:** 36 total items

**Models:** `gpt-5-nano` first-pass, `gpt-5.2` re-rank top 40

---

## Neuroblastoma (21 shown / 21 total)

### [Whole genome sequencing of endometrial cancer identifies novel subgroups, drivers, and actionable alterations](https://www.biorxiv.org/content/10.1101/2026.06.15.730391v1?rss=1)
*bioRxiv Genomics*  
Score: **0.85**
Published: 2026-06-19
Tags: endometrial cancer, WGS, biomarkers, actionable

Preprint reporting endometrial cancer subgroups and actionable genomic alterations, with clear translational biomarker potential.

---

### [Accurate detection of tumor clonality and ongoing expansion mode from genomic data](https://www.biorxiv.org/content/10.64898/2026.06.15.732415v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.66**
Published: 2026-06-19
Tags: tumor evolution, clonality, intra-tumor heterogeneity, mutation clustering, resistance

Introduces DECODE, a mutation clustering method that models sample-specific coverage and calling biases to improve intra-tumor heterogeneity/clonality estimates. Useful for studying clonal dynamics that often underpin treatment resistance and relapse, though cancer-type context is not specified here.

<details>
<summary>RSS summary</summary>

Recent evidence shows that despite considerable effort, currently available algorithms for estimating intra-tumor heterogeneity (ITH) remain limited. We developed DECODE (Deciphering Cancer Origin from DNA Evolution), a novel mutation clustering method that incorporates the impact of sample-specific sequencing coverage and mutation calling biases. On synthetic data, DECODE outperformed existing methods across multiple clonality metrics and accurately detected and characterized the neutral tail i…

</details>

---

### [From hotspot dependence to distributed robustness in resistance-aware lead optimization](https://www.biorxiv.org/content/10.64898/2026.06.16.732538v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.55**
Published: 2026-06-22
Tags: drug resistance, lead optimization, targeted therapy, mutational liabilities, computational drug design

Proposes a framework (ResistAgent) to incorporate mutational resistance liabilities into lead optimization with resistance mapping and counter-design. Directly about resistance-aware drug design, but examples are EGFR/HIV and not pediatric oncology or omics-focused.

<details>
<summary>RSS summary</summary>

Drug resistance remains a recurrent failure mode in targeted anticancer and antiviral therapy, and resistance evidence often enters only after compound selection. ResistAgent is an evidence-constrained framework that converts mutational liabilities into design-time objectives through site- and combo-aware resistance mapping, deterministic mechanism diagnosis and robust counter-design. In EGFR-Erlotinib and HIV-RT-Rilpivirine, the framework separated residue-level liabilities from observed HIV co…

</details>

---

### [Argonaute 2 drives resistance to immune checkpoint inhibitors in immunorefractory non-small cell lung cancer](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3003860)
*PLOS Biol*  
Score: **0.46**
Published: 2026-06-18
Tags: therapy resistance, immunotherapy, immune checkpoint inhibitors, NSCLC, mechanism

Mechanistic study of immune checkpoint inhibitor resistance in NSCLC implicating Argonaute 2, which is relevant to therapy resistance biology. However, it is not pediatric-focused and the summary does not indicate a computational or multi-omics methods angle.

<details>
<summary>RSS summary</summary>

<p>by Dario Pasquale Anobile, Layla Barbar, Emile Maucotel, Alexis Cornec, Valeria Manriquez, Wilfrid Richer, Jordan Denizeau, Christine Sedlik, Charlie Bories, Elodie Couderc, Renaud Leclere, Judith Sobas, Emeline Papillon, Rafael Mena Osuna, Jimena Tosello-Boari, Marianne Burbage, Eliane Piaggio, Enzo Z. Poirier</p> One of the first-line treatments for advanced non-small cell lung cancer (NSCLC) are immune checkpoint inhibitors (ICI), which activate the antitumor immune response. Despite their…

</details>

---

### [Additivity, Not Synergy, Underlies the Efficacy of Current Combination Regimens in Urothelial Cancer](https://aacrjournals.org/cancerrescommun/article/6/6/1447/785841/Additivity-Not-Synergy-Underlies-the-Efficacy-of)
*AACR Cancer Research Communications*  
Score: **0.45**
Published: 2026-06-19
Tags: combination therapy, additivity vs synergy, immunotherapy, clinical trials, translational

Analyzes whether ICI-based combinations are synergistic vs additive in urothelial cancer, a translationally relevant framing for combination strategy evaluation. Not pediatric/NB and the summary does not indicate omics-based biomarkers or single-cell methodology.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle"></div>PD-1/PD-L1 immune checkpoint inhibitor (ICI)-based combination regimens have transformed metastatic urothelial cancer treatment. However, phase III trials of these regimens have shown varied results. We sought to understand these differences by assessing whether these combinations are synergistic, additive, or less-than-additive, which is crucial for translational research and clinical development. We …

</details>

---

### [OracleScreen-LILRB4: Machine Learning-Guided Discovery of Myeloid Immune Checkpoint Binders Validated in Patient-Derived Cells](https://www.biorxiv.org/content/10.64898/2026.06.17.732859v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.41**
Published: 2026-06-21
Tags: machine learning, drug discovery, immune checkpoint, LILRB4, patient-derived cells

Uses an ensemble ML framework to identify small-molecule modulators of an immune checkpoint (LILRB4) with validation in patient-derived cells. Translational drug-discovery oriented, but not tied to pediatric oncology or single-cell/multi-omics analysis workflows.

<details>
<summary>RSS summary</summary>

The identification of small molecule modulators of immune checkpoint proteins remains a significant challenge in drug discovery due to the flat, featureless nature of protein-protein interaction interfaces and the characteristically low hit rates observed in conventional high-throughput screening campaigns. Here we report OracleScreen-LILRB4, an ensemble machine learning framework trained on quantitative biophysical screening data from two structurally diverse compound libraries (19,800 compound…

</details>

---

### [Cancer cells adopt unprecedented strategies to produce a molecule that protects them from iron-dependent death](https://www.nature.com/articles/d41586-026-01802-3)
*Nature*  
Score: **0.41**
Published: 2026-06-22
Tags: cancer, ferroptosis, iron metabolism, cell death, therapy resistance

Discusses spermine binding iron to prevent ferroptosis, a cancer-relevant cell-death/resistance concept; however this appears to be a news/commentary-style piece with limited computational or pediatric translational detail in the summary.

<details>
<summary>RSS summary</summary>

<p>Nature, Published online: 22 June 2026; <a href="https://www.nature.com/articles/d41586-026-01802-3">doi:10.1038/d41586-026-01802-3</a></p>The finding that spermine molecules in cells bind to iron to prevent it unleashing ferroptosis, a type of cell death, opens up strategies for treating tissue damage and cancer.

</details>

---

### [HTS-Oracle v2: Prospective AI-Guided Discovery and Experimental Validation of Small Molecule Modulators Across Multiple Targets](https://www.biorxiv.org/content/10.64898/2026.06.15.732399v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.38**
Published: 2026-06-19
Tags: AI, high-throughput screening, drug discovery, immune checkpoints, model validation

AI-guided hit discovery for immune checkpoint targets in HTS with emphasis on rigorous cross-validation and prospective validation. Useful conceptually for drug discovery, but limited direct relevance to single-cell/multi-omics biomarker discovery in pediatric oncology.

<details>
<summary>RSS summary</summary>

High-throughput screening (HTS) remains the cornerstone of early-phase small molecule discovery yet consistently underperforms against immunotherapy targets, yielding validated hit rates below 0.1%. Here we introduce HTS-Oracle v2, which features rigorous cross-validation that ensures honest performance estimates. HTS-Oracle v2 was trained and validated across four clinically significant immune checkpoint targets (CD28, ICOS, LAG-3, and TIGIT) achieving ROC-AUC values of 0.968, 0.969, 0.875, 0.9…

</details>

---

### [Beyond the canonical: The role of post-transcriptional regulation in drug-target interaction prediction](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014440)
*PLOS Comp Bio*  
Score: **0.33**
Published: 2026-06-22
Tags: drug-target interaction, alternative splicing, isoforms, computational pharmacology

Highlights how alternative splicing/isoforms can affect drug binding and proposes incorporating post-transcriptional regulation into DTI/DTA prediction. Relevant to drug response modeling in principle, but not clearly connected to cancer cohorts, biomarkers, or single-cell omics in the provided summary.

<details>
<summary>RSS summary</summary>

<p>by Md Istiaq Ansari, Khandakar Tanvir Ahmed, Debby D. Wang, Kirill Medvedev, Wei Zhang</p> Protein isoforms produced from the same gene through post-transcriptional regulatory mechanisms, such as alternative splicing, can substantially alter protein structure and function, including drug-binding properties. However, most existing drug-target interaction (DTI) and drug-target affinity (DTA) prediction models rely exclusively on a single representative protein sequence per gene, typically the c…

</details>

---

### [Geometric Deep Learning Reveals Ligandable and Cryptic RNA Binding Small Molecule Pockets (SMARTPocket)](https://www.biorxiv.org/content/10.64898/2026.06.18.732920v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.28**
Published: 2026-06-19
Tags: geometric deep learning, RNA, drug discovery, binding pockets, transfer learning

Geometric deep learning method to predict RNA small-molecule binding pockets from 3D structure with transfer learning from protein interfaces. Primarily RNA-targeted drug discovery and not tied to pediatric oncology, biomarker discovery, or single-cell omics.

<details>
<summary>RSS summary</summary>

RNAs are important therapeutic targets, however identifying ligandable small-molecule binding pockets remains a major barrier to RNA-targeted drug discovery. Here, SMARTPocket, an atomic-level geometric deep learning framework for predicting RNA-small molecule binding pockets directly from three-dimensional structure is introduced. SMARTPocket represents RNA as full-atom point clouds and uses transfer learning from more than 110,000 protein binding interface structures to overcome the limited nu…

</details>

---

### [Integrative modelling of innate immune response dynamics during virus infection](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014395)
*PLOS Comp Bio*  
Score: **0.24**
Published: 2026-06-22
Tags: systems biology, dynamic modeling, innate immunity, viral infection, integrative modeling

Quantitative/integrative dynamic modeling of host–virus innate immune interactions may be broadly informative for systems modeling, but it is not tied to cancer, biomarkers, or single-cell omics in the provided summary.

<details>
<summary>RSS summary</summary>

<p>by Ramya Boddepalli, Harsh Chhajera, Rahul Roya</p> Positive-sense RNA viruses that constitute a large class of human pathogens employ various strategies to suppress and evade host immune defenses. Understanding the dynamic interaction between the viral life cycle and immune signaling is crucial to designing effective antiviral strategies. Although significant progress has been made, quantitative models that can accurately capture the intricate interactions and the intertwined dynamics during…

</details>

---

### [CD4+ T cells impair tumor growth through IL-3 and TNF-dependent vascular damage](https://www.science.org/doi/abs/10.1126/science.ads7910?af=R)
*Science*  
Score: **0.21**
Published: 2026-06-18
Tags: tumor immunology, CD4 T cells, cytokines, vascular

Immunology/mechanism paper on CD4+ T cells affecting tumor growth via cytokines and vascular damage. Potential background for tumor-immune interactions, but no single-cell/omics or pediatric/translational biomarker emphasis is evident from the provided entry.

<details>
<summary>RSS summary</summary>

Science, Volume 392, Issue 6804, June 2026. <br />

</details>

---

### [HoloBio: A holographic microscopy tool for quantitative biological analysis](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013928)
*PLOS Comp Bio*  
Score: **0.10**
Published: 2026-06-22
Tags: microscopy, image analysis, software, holographic imaging

Software for quantitative analysis of holographic microscopy data. General bioimage analysis utility, but not clearly connected to oncology, omics, or biomarker/resistance themes.

<details>
<summary>RSS summary</summary>

<p>by Waira Mona, Maria J. Gil-Herrera, Emanuel Mazo, Daniel Córdoba, Sofia Obando-Vasquez, Maria J. Lopera, Rene Restrepo, Carlos Trujillo, Ana Doblas, Raul Castaneda</p> Holographic imaging in microscopy enables label-free quantitative information of biological specimens and has found applications across a wide range of biomedical studies, from cell morphology to particle dynamics; yet its widespread adoption is often limited by the lack of accessible and standardized analysis software. We pre…

</details>

---

### [GrassSV – hybrid method to detect structural variants in high throughput DNA-seq data](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014406)
*PLOS Comp Bio*  
Score: **0.09**
Published: 2026-06-22
Tags: structural variants, DNA-seq, variant calling, method

Structural-variant detection method framed around genetic diversity in populations, with no explicit cancer or clinical genomics angle in the summary. Could be tangentially useful, but relevance to translational pediatric oncology is unclear from the provided text.

<details>
<summary>RSS summary</summary>

<p>by Dominik Witczak, Krzysztof Sychla, Julia Wysocka, Artur Laskowski, Wojciech Frohmberg, Marta Glowacka, Alicja Dzik, Piotr Lukasiak, Jacek Blazewicz, Aleksandra Swiercz</p> Genetic diversity is crucial for populations to adapt and survive in dynamic environments. This diversity arises from genetic mutations, which manifest in the genome as structural variants (SVs). Several types of SVs exist, but not all are equally easy to detect. Current SV detection tools tend to specialize in certain S…

</details>

---

### [A comparison of contact patterns derived from the population structure in agent-based models and empirical contact survey data](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013533)
*PLOS Comp Bio*  
Score: **0.08**
Published: 2026-06-18
Tags: agent-based model, contact patterns, infectious disease modeling

Focuses on agent-based modeling of infectious-disease contact patterns and comparison to survey data, which is outside pediatric oncology and omics/single-cell method priorities.

<details>
<summary>RSS summary</summary>

<p>by Janik Suer, Johannes Ponge, Michael Brüggemann, Jan Pablo Burgard, Vitaly Belik, Bernd Hellingrath, Alejandra Rincón Hidalgo, Andrzej K. Jarynowski, Richard Pastor, Huynh Thi Phuong, Steven Schulz, Ashish Thampi, Chao Xu, Marlli Zambrano, Rafael Mikolajczyk, André Karch, Veronika K. Jaeger, on behalf of the OptimAgent Consortium </p> Agent-based models (ABMs) are powerful tools for simulating disease spread, relying on individual-level interaction rules from which emergent dynamics arise. …

</details>

---

### [Big Ebola outbreak puts research spotlight on little-known virus](https://www.science.org/doi/abs/10.1126/science.aej7944?af=R)
*Science*  
Score: **0.08**
Published: 2026-06-18
Tags: virology, public health

Virology outbreak coverage; limited translational relevance to NB.

---

### [Mechanisms underlying spontaneous and evoked calcium responses in oligodendrocyte precursor cells: A modeling investigation](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013430)
*PLOS Comp Bio*  
Score: **0.07**
Published: 2026-06-18
Tags: calcium signaling, computational modeling, oligodendrocytes, neuroscience

Computational modeling of Ca2+ signaling in oligodendrocyte precursor cells is neuroscience-focused and not connected to cancer, biomarkers, therapy response, or omics integration in the summary.

<details>
<summary>RSS summary</summary>

<p>by Martin Lardy, Leqi Wang, Claire Guerrier, Veronica T. Cheli, Pablo M. Paez, Anmar Khadra</p> Calcium (Ca<sup>2+</sup>) signaling has emerged as a central regulator of activity-dependent myelination in oligodendrocytes. These Ca<sup>2+</sup> signals encompass both the stimulus-independent spontaneous Ca<sup>2+</sup> local transients (SCaLTs) generated intrinsically in a voltage-independent manner or facilitated by the membrane voltage, as well as evoked responses triggered by ATP and glutam…

</details>

---

### [Developmental oligodendrocytes regulate brain function through the mediation of synchronized spontaneous activity](https://elifesciences.org/articles/102200)
*eLife*  
Score: **0.05**
Published: 2026-06-18
Tags: neuroscience, oligodendrocytes, development

Developmental neuroscience study on oligodendrocytes and synchronized activity with no oncology, biomarker, or computational multi-omics emphasis evident from the summary snippet.

<details>
<summary>RSS summary</summary>

Synchronized spontaneous neural activity is a fundamental feature of developing central nervous systems and is thought to be essential for proper brain development. However, the mechanisms that regulate this synchronization and its long-term impact on brain function remain unclear. Here, we identify a previously unrecognized role of oligodendrocytes in orchestrating synchronized spontaneous activity during a critical developmental window, with lasting consequences for adult behavior. Using oligo…

</details>

---

### [Enterovirus D68 2A protease causes nuclear pore complex dysfunction and independently contributes to motor neuron toxicity](https://elifesciences.org/articles/108672)
*eLife*  
Score: **0.05**
Published: 2026-06-18
Tags: virology, motor neuron, nuclear pore complex, EV-D68

Mechanistic virology/neurotoxicity work on EV-D68 and nuclear pore complex dysfunction; not aligned with pediatric oncology, resistance, or omics method development based on the summary.

<details>
<summary>RSS summary</summary>

Enterovirus D68 (EV-D68) is an important pathogen associated with acute flaccid myelitis (AFM). The pathogenesis of AFM involves infection of spinal motor neurons and motor neuron death; however, the mechanisms linking EV-D68 infection to selective neurotoxicity are not well understood. Dysfunction of the nuclear pore complex (NPC) has been implicated in motor neuron injury in neurodegenerative diseases such as amyotrophic lateral sclerosis, and the NPC is also modified by picornavirus proteases…

</details>

---

### [Induction of broadly neutralizing HIV antibodies by a two-step mechanism informs vaccine design](https://www.science.org/doi/abs/10.1126/science.aec6396?af=R)
*Science*  
Score: **0.02**
Published: 2026-06-18
Tags: HIV, vaccine design, neutralizing antibodies

HIV vaccine/antibody induction study; not aligned with pediatric oncology, resistance, or biomarker/multi-omics priorities based on the citation-only summary.

<details>
<summary>RSS summary</summary>

Science, Volume 392, Issue 6804, June 2026. <br />

</details>

---

### [Artificial light pollution disrupts sleep and neuronal genomic stability in wild reef fish](https://www.cell.com/current-biology/fulltext/S0960-9822(26)00663-9?rss=yes)
*Cell*  
Score: **0.00**
Published: 2026-06-22
Tags: ecology, sleep, DNA damage, reef fish

Ecology/neurobehavior study in reef fish about light pollution and neuronal DNA damage; outside cancer, pediatric oncology, and omics computational method scope.

<details>
<summary>RSS summary</summary>

Ben-Ezra et al. use lab and field experiments to show that artificial light at night alters territorial occupancy, aggression, feeding, and sleep in Chromis viridis. Acute and chronic light pollution-driven sleep loss correlates with increased neuronal DNA damage in the dorsal pallium of these reef fish.

</details>

---

## AI (2 shown / 2 total)

### [EventHorizon: A Foundation Model for Clinical Flow Cytometry](https://www.biorxiv.org/content/10.64898/2026.06.18.733197v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.84**
Published: 2026-06-22
Tags: foundation model, flow cytometry, self-supervised, representation learning, clinical translational, biomarkers

Presents a self-supervised foundation model to learn unified specimen-level representations from clinical flow cytometry, explicitly addressing generalization across panels/instruments/labs. Strongly aligned with representation learning for clinically deployed single-cell(-like) assays and potential biomarker/diagnostic workflows.

<details>
<summary>RSS summary</summary>

Flow cytometry is an essential tool for diagnosis of hematologic malignancies, but existing clinical workflows are highly dependent on expert manual interpretation. Existing machine learning approaches typically require extensive labeled data and are sensitive to variability in panel design, instrumentation, and laboratory workflows, limiting their generalizability. We present EventHorizon, a self-supervised foundation model for clinical flow cytometry that produces unified specimen-level repres…

</details>

---

### [Morpho-FM: spatial molecular reconstruction from routine H&E histology using transcriptomic foundation-model priors](https://www.biorxiv.org/content/10.64898/2026.06.15.732498v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.73**
Published: 2026-06-19
Tags: computational pathology, H&E, spatial reconstruction, foundation model, transcriptomics, translational

Uses transcriptomic foundation-model priors to infer spatial gene-expression programs from routine H&E, addressing small paired-cohort limitations of prior histology-to-expression models. Potentially impactful for translational pathology workflows where spatial assays are sparse/costly.

<details>
<summary>RSS summary</summary>

Routine haematoxylin and eosin (H&amp;E) histology captures tissue architecture at clinical scale, but lacks a direct molecular readout of the transcriptional programmes that organise tumour epithelium, stroma, vasculature and immune compartments. Spatial transcriptomics provides this context, yet cost, workflow complexity and sparse sampling limit routine use. Most existing histology-to-expression models are trained de novo on small paired cohorts and therefore remain weakly constrained when ex…

</details>

---

## Methods (2 shown / 2 total)

### [scMagnifier: Resolving fine-grained cell subtypes via GRN-informed perturbations and consensus clustering](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014167)
*PLOS Comp Bio*  
Score: **0.78**
Published: 2026-06-18
Tags: scRNA-seq, cell states, clustering, gene regulatory network, regulatory inference, heterogeneity

Introduces a scRNA-seq consensus clustering method that uses gene regulatory network–informed in silico perturbations to amplify subtle expression differences and reveal latent subpopulations. Useful for resolving tumor cell-state heterogeneity and regulatory programs tied to response/resistance.

<details>
<summary>RSS summary</summary>

<p>by Zhenhui He, Dong Kangning</p> Resolving fine-grained cell subtypes in single-cell RNA sequencing (scRNA-seq) data remains challenging, as their subtle transcriptional differences are often obscured by technical noise and data sparsity. Here, we present scMagnifier, a consensus clustering framework that leverages gene regulatory network (GRN)-informed in silico perturbations to amplify subtle transcriptional differences and uncover latent cell subpopulations. scMagnifier perturbs candidate …

</details>

---

### [Lamprey 3D single-cell transcriptomics reveals ancestral and specialized features of the vertebrate brain](https://www.science.org/doi/abs/10.1126/science.aea2535?af=R)
*Science*  
Score: **0.30**
Published: 2026-06-18
Tags: single-cell transcriptomics, spatial/3D, brain, methods-adjacent

Uses 3D single-cell transcriptomics, which is methodologically adjacent to scRNA-seq workflows, but the biological focus is vertebrate brain evolution rather than oncology or translational biomarker discovery (no further details in the summary).

<details>
<summary>RSS summary</summary>

Science, Volume 392, Issue 6804, June 2026. <br />

</details>

---

## Other (11 shown / 11 total)

### [SPOmiAlign: a modality-agnostic computational framework for multimodal spatial omics alignment enabled by a feature matching foundation model](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag331/8712668?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.82**
Published: 2026-06-21
Tags: spatial omics, multi-omics, alignment, foundation model, data integration, translational

Targets a core bottleneck in spatial multi-omics: robust alignment across modalities/sections, using a feature-matching foundation-model approach. Directly relevant to integrating transcriptomic/proteomic/metabolomic spatial readouts for tumor microenvironment and state mapping.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Multimodal spatial omics enables systematic characterization of tissue organization by jointly profiling transcriptomic, proteomic, metabolomic, and other spatially resolved modalities within their spatial context. A central challenge in realizing this potential is achieving robust spatial alignment across modalities and sections. Although numerous alignment methods have been developed, most are designed for single-modality secti…

</details>

---

### [ModelistsGCN: a multimodal graph convolutional network framework for single-cell spatial transcriptomic cell typing](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag340/8713041?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.74**
Published: 2026-06-22
Tags: spatial transcriptomics, cell typing, graph neural network, multimodal, single-cell

Proposes a multimodal GCN to improve cell typing in single-cell spatial transcriptomics where gene panels are limited. Directly relevant to spatially resolved tumor ecosystems and translating scRNA-seq labels to spatial platforms.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>Spatial transcriptomics technologies currently face a trade-off between spatial and molecular resolution. High-resolution single-cell methods such as MERFISH and Expansion Sequencing (ExSeq) resolve individual cells but typically profile limited gene panels, restricting the direct application of transcriptome-based cell typing strategies developed for single-cell RNA sequencing. Consequently, accurate cell-type identification in …

</details>

---

### [ST2HE: enhancing spatial transcriptomics interpretability via virtual staining for histological annotation](https://academic.oup.com/bib/article/doi/10.1093/bib/bbag326/8713038?rss=1)
*Briefings in Bioinformatics (Oxford Academic)*  
Score: **0.67**
Published: 2026-06-22
Tags: spatial transcriptomics, virtual staining, diffusion model, interpretability, histology annotation

Presents a diffusion-based generative method to synthesize virtual H&E images from high-resolution spatial transcriptomics to aid histological annotation. Relevant for improving interpretability and linking molecular states to histology in translational tumor studies.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div>High-resolution spatial transcriptomics (HR-ST) technologies offer unprecedented insights into tissue architecture but lack standardized frameworks for histological annotation. We present ST2HE, a cross-platform generative framework that synthesizes virtual hematoxylin and eosin images directly from HR-ST data. ST2HE integrates nuclei morphology and spatial transcript coordinates using a one-step diffusion model, enabling histolo…

</details>

---

### [ChromBERT-tools: A versatile toolkit for context-specific regulatory representations of transcription regulators across different cell types](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btag423/8711926?rss=1)
*Bioinformatics (Oxford Academic)*  
Score: **0.63**
Published: 2026-06-19
Tags: regulatory representation, transcription regulators, BERT, gene regulation, toolkit

Provides a toolkit for learning context-specific regulatory representations of transcription regulators across cell types, aiming to support flexible transcription modeling and in silico regulatory analysis. Potentially transferable to regulatory-program inference for tumor cell states, but the summary is not cancer-specific.

<details>
<summary>RSS summary</summary>

<span class="paragraphSection"><div class="boxTitle">Abstract</div><div class="boxTitle">Summary</div>Representations that encode the genome-wide regulatory behavior of transcription regulators provide a foundation for flexible transcription modeling and in silico regulatory analysis. Existing regulator representations are commonly derived from gene co-expression, motif annotations, or static protein features, which capture useful but limited aspects of regulator identity but do not directly mod…

</details>

---

### [Antibody-Antigen Affinity Prediction with Chain-Aware Protein Language Modeling](https://www.biorxiv.org/content/10.64898/2026.06.19.733375v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.52**
Published: 2026-06-21
Tags: protein language model, representation learning, antibody, affinity prediction, bioinformatics, transfer learning

Introduces a sequence-based, chain-aware protein language modeling approach for predicting antibody–antigen affinity, which is a transferable representation-learning method even though it is not cancer or single-cell focused in the provided summary.

<details>
<summary>RSS summary</summary>

Motivation: Antibody-antigen affinity determines which antibodies advance in therapeutic discovery, repertoire analysis and affinity maturation, but experimental measurements are sparse relative to the scale of sequence libraries. Structure-based predictors can exploit interface geometry when reliable complexes are available, yet early discovery often requires ranking many heavy-light chain pairs against antigens for which no complex structure exists. Existing sequence-based models are scalable,…

</details>

---

### [CoDaLoMic: An R package for modeling microbiome compositional and longitudinal data](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014328)
*PLOS Comp Bio*  
Score: **0.43**
Published: 2026-06-22
Tags: longitudinal modeling, compositional data, R package, statistical methods, omics methods

Presents an R package with models tailored to longitudinal, compositional data and claims to capture interactions among groups of bacteria, offering statistical tooling that could transfer to other longitudinal omics settings despite being microbiome-focused.

<details>
<summary>RSS summary</summary>

<p>by Irene Creus-Martí, Andrés Moya, Francisco J. Santonja</p> In this paper we present <i>CoDaLoMic</i>, an R package for analyzing longitudinal and compositional microbiome datasets. The <i>CoDaLoMic</i> package implements three models specifically designed for the analysis of microbiome data that are both compositional and longitudinal. Unlike many existing methods that focus solely on pairwise interactions, <i>CoDaLoMic</i> also captures interactions among groups of bacteria, providing a mo…

</details>

---

### [TCRBinder: Unified pre-trained language model with paired-chain synergy for predicting T-cell receptor binding specificity](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014396)
*PLOS Comp Bio*  
Score: **0.37**
Published: 2026-06-22
Tags: protein language model, TCR, pHLA, immunology, pretraining

Develops a pre-trained language model for predicting TCR specificity using paired-chain information. Methodologically strong but more immunology/antigen-specificity focused than tumor multi-omics/single-cell state modeling described in the interests.

<details>
<summary>RSS summary</summary>

<p>by Weihe Dong, Qiang Yang, Long Xu, Xiaokun Li, Kuanquan Wang, Suyu Dong, Gongning Luo, Xianyu Zhang, Tiansong Yang, Xin Gao, Guohua Wang</p> Deciphering how human T cells recognise peptide-HLA (pHLA) complexes underpins next-generation vaccines and personalised immunotherapies, yet extreme sequence diversity and paired-chains interdependence still hamper reliable <i>in silico</i> prediction of T-cell receptor (TCR) specificity. To overcome these hurdles, we built TCRBinder, a paired-chain-aw…

</details>

---

### [SAbDab2: The structural antibody database in the age of machine learning](https://www.biorxiv.org/content/10.64898/2026.06.16.732554v1?rss=1)
*bioRxiv Bioinformatics*  
Score: **0.22**
Published: 2026-06-20
Tags: database, antibodies, structures, machine learning resource

Announces an updated structural antibody database designed for modern ML use and new antibody formats. Useful infrastructure for antibody modeling, but distant from single-cell/multi-omics translational oncology workflows in the prompt.

<details>
<summary>RSS summary</summary>

The Structural Antibody Database (SAbDab) is a publicly available repository of experimentally determined antibody structures, first released in 2013. Explicit support for single-domain antibodies was added in 2021, with SAbDab-nano. Recently, increasing interest in antibodies has led to a proliferation of novel antibody formats, while simultaneous advances in machine learning have increased demand for standardised, high-quality structure data. Here, we present SAbDab2, re-engineered for the mac…

</details>

---

### [Maternal trans-vaccenic acid shapes neonatal T cell development and early-life immune imprinting](https://www.science.org/doi/abs/10.1126/science.aea4041?af=R)
*Science*  
Score: **0.04**
Published: 2026-06-18
Tags: immunology, T cells, development, metabolite

Immunology/developmental study on neonatal T cell development; lacks cancer, translational oncology, or computational omics methods in the provided information.

<details>
<summary>RSS summary</summary>

Science, Volume 392, Issue 6804, June 2026. <br />

</details>

---

### [Neural signatures of model-based and model-free reinforcement learning across prefrontal cortex and striatum](https://elifesciences.org/articles/106032)
*eLife*  
Score: **0.03**
Published: 2026-06-22
Tags: neuroscience, reinforcement learning, prefrontal cortex, striatum

Reinforcement learning/neurophysiology paper focused on PFC/striatum; not relevant to translational oncology or omics/single-cell computational biology from the summary.

<details>
<summary>RSS summary</summary>

Animals integrate knowledge about how the state of the environment evolves to choose actions that maximise reward. Such goal-directed behaviour – or model-based (MB) reinforcement learning (RL) – can flexibly adapt choice to changes, being thus distinct from simpler habitual – or model-free (MF) RL – strategies. Previous inactivation and neuroimaging work implicates prefrontal cortex (PFC) and the caudate striatal region in MB-RL; however, details are scarce about its implementation at the singl…

</details>

---

### [Brain-wide topographic coordination of rotating waves](https://www.science.org/doi/abs/10.1126/science.adx1369?af=R)
*Science*  
Score: **0.02**
Published: 2026-06-18
Tags: neuroscience, brain dynamics

Systems neuroscience topic (brain-wide rotating waves) with no apparent connection to oncology or computational omics methods in the citation-only summary.

<details>
<summary>RSS summary</summary>

Science, Volume 392, Issue 6804, June 2026. <br />

</details>

---
