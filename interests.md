# tocify
# Ingest ToCs and push links to papers relevant to me

# Interests seed

## Keywords
EEG
MEG
neural oscillations
aperiodic
neural timescales
ECG
respiration
parameterization
specparam
fooof
bycycle

## Current narrative
# These are just from my NIH biosketch
Neural data science: My research program makes extensive use of many heterogeneous datasets
from multiple different species in order to converge on our physiological and cognitive questions of
interest. We leverage large-scale data analysis and we have moved to an “open science” format,
including the publishing of available data (when appropriate and authorized) as well as publishing all
analysis code using free and open software (e.g., Python). I’ve demonstrated how data mining of the
peer-reviewed neuroscientific literature can aid research and scientific discovery (b). In an invited
commentary, I argue the advantages of creating an open data ecosystem for facilitating these kinds of
data-mining projects (a) and (c). All of my methods and data (where possible) are made public, and
my lab has been developing interactive tutorials and publishing software packages made freely
available to everyone on GitHub, including software for basic neural data digital signal processing. I
have adopted an open education approach in my undergraduate data science teaching, leveraging a
variety of technological platforms, which we have recently formalized and shared (d). In an upcoming,
accepted Commentary to be published in Nature Methods, I argue that the future of Neural Data Science will be in combing the many, massive, heterogeneous datasets being collected around the
world (human fMRI, gene expression, connectivity, etc.). And that, once integrated, these datasets
can be mined across features, which will help spur the development of new hypotheses.
a. Donoghue T, Voytek B, Ellis S. Teaching Creative and Practical Data Science at Scale. Journal of
Statistics and Data Science Education. 2021 March 22; 29(sup1):S27-S39.
b. Voytek B. Social Media, Open Science, and Data Science Are Inextricably Linked. Neuron. 2017
Dec 20;96(6):1219-1222. PubMed PMID: 29224728.
c. Voytek B. The Virtuous Cycle of a Data Ecosystem. PLoS Comput Biol. 2016
Aug;12(8):e1005037. PubMed Central PMCID: PMC4974004.
d. Voytek JB, Voytek B. Automated cognome construction and semi-automated hypothesis
generation. J Neurosci Methods. 2012 Jun 30;208(1):92-100. PubMed Central PMCID:
PMC3376233.

Neural oscillation waveform shape informs physiology and cognition: My lab has recently pioneered
an entirely new physiological framework for analyzing the shape of neural oscillations, showing that
nonsinusoidal features of oscillations may reflect alterations in synaptic input in e.g., Parkinson’s
disease (b). Additionally, we introduced a novel suite of tools for analyzing the shape of local field
potential and EEG waveforms, which we argue provide for better inference of the physiology of the
field potential generators and how they are influenced by cognitive, behavioral, and disease states (a).
Next, we have introduced a suite of Python tools—tested and validated against simulated ground-truth
and animal data—for performing these analyses, including full tutorials, demonstrating how this
approach gives better estimates of oscillation features (c). Leveraging our waveform shape analysis
and spectral parameterization toolboxes funded by our current NIGMS R01, and a longitudinal dataset
of publicly-available baby EEG data, we showed that oscillation waveform shape and aperiodic activity
develops quickly in infants in the first few months after birth (d).
a. Schaworonkow N, Voytek B. Longitudinal changes in aperiodic and periodic activity in
electrophysiological recordings in the first seven months of life. Dev Cogn Neurosci. 2021
Feb;47:100895. PubMed Central PMCID: PMC7734223.
b. Cole S, Voytek B. Cycle-by-cycle analysis of neural oscillations. J Neurophysiol. 2019 Aug
1;122(2):849-861. PubMed PMID: 31268801.
c. Cole SR, van der Meij R, Peterson EJ, de Hemptinne C, Starr PA, Voytek B. Nonsinusoidal Beta
Oscillations Reflect Cortical Pathophysiology in Parkinson's Disease. J Neurosci. 2017 May
3;37(18):4830-4840. PubMed Central PMCID: PMC5426572.
d. Cole SR, Voytek B. Brain Oscillations and the Importance of Waveform Shape. Trends Cogn Sci.
2017 Feb;21(2):137-149. PubMed PMID: 28063662.

Aperiodic activity tracks excitation/inhibition (EI) balance and cognition: We previously showed that
1/f-like spectral “neural noise” is associated with cognitive impairments and age-related working
memory decline. More recently we have explored the neural origins of changes in this 1/f-like aspect
of the power spectrum. Specifically, we developed a novel theoretical framework whereby we argued
that EI imbalances impact neural spiking and information transmission (b). Next, using a mix of
computational modeling as well as data from rat hippocampal recordings and macaque cortex, we
demonstrated that this aperiodic, 1/f-like signal may index EI balance (a). Most recently, with our
current NIGMS R01 support, we developed an extensive analytic framework for disentangling this
aperiodic activity from oscillations, providing a robust parameterization toolbox, demonstrating how it
can be leveraged to uncover novel links to working memory, and to analyze very large datasets (d).
Leveraging this toolbox and a longitudinal dataset of baby EEG, we showed that aperiodic activity
develops quickly in infants in the first few months after birth (c).
a. Schaworonkow N, Voytek B. Longitudinal changes in aperiodic and periodic activity in
electrophysiological recordings in the first seven months of life. Dev Cogn Neurosci. 2021
Feb;47:100895. PubMed Central PMCID: PMC7734223.
b. Donoghue T, Haller M, Peterson EJ, Varma P, Sebastian P, Gao R, Noto T, Lara AH, Wallis JD,
Knight RT, Shestyuk A, Voytek B. Parameterizing neural power spectra into periodic and aperiodic
components. Nat Neurosci. 2020 Dec;23(12):1655-1665. PubMed Central PMCID: PMC8106550.
c. Gao R, Peterson EJ, Voytek B. Inferring synaptic excitation/inhibition balance from field potentials.
Neuroimage. 2017 Sep;158:70-78. PubMed PMID: 28676297.
d. Voytek B, Knight RT. Dynamic network communication as a unifying neural basis for cognition,
development, aging, and disease. Biol Psychiatry. 2015 Jun 15;77(12):1089-97. PubMed Central
PMCID: PMC4443259.

Oscillations as neural computational mechanism: To examine the role that neural oscillations play in
cognition and in coordinating information transfer between brain regions, we have made significant
use of invasive electroencephalography (iEEG) in people undergoing surgery for intractable epilepsy.
To address our cognitive questions of interest, we have developed several novel methods, including a
novel metric for assessing the relationship between low-frequency oscillations and local population
spiking activity, known as phase-amplitude coupling (PAC) (b). Taking advantage of this new
approach, we showed how oscillatory dynamics play a critical role in coordinating cognitive networks
during human goal maintenance (a). We found that oscillatory networks form quickly, with different
task goals represented at different coupling phases to permit directional information transfer on a trial-
by-trial basis. Through collaboration, we have examined the physiological origin of oscillations and
their development through the use of 3D, human induced pluripotent stem cell cortical organoids. We
found that oscillations naturally emerge with development, dependent on functioning excitatory /
inhibitory (glutamatergic and GABAergic) interactions (c). Finally, to understand the potential
functional / computational role that oscillations play in cognition, we must first be certain we are
actually measuring oscillations as intended. To that end, with our current NIGMS R01 support, we
recently published a checklist of seven technical points that one must carefully consider when
approaching the study of oscillations (d).
a. Donoghue T, Schaworonkow N, Voytek B. Methodological considerations for studying neural
oscillations. Eur J Neurosci. 2022 Jun;55(11-12):3502-3527. PubMed Central PMCID:
PMC8761223.
b. Trujillo CA, Gao R, Negraes PD, Gu J, Buchanan J, Preissl S, Wang A, Wu W, Haddad GG,
Chaim IA, Domissy A, Vandenberghe M, Devor A, Yeo GW, Voytek B, Muotri AR. Complex
Oscillatory Waves Emerging from Cortical Organoids Model Early Human Brain Network
Development. Cell Stem Cell. 2019 Oct 3;25(4):558-569.e7. PubMed Central PMCID:
PMC6778040.
c. Voytek B, Kayser AS, Badre D, Fegen D, Chang EF, Crone NE, Parvizi J, Knight RT, D'Esposito
M. Oscillatory dynamics coordinating human frontal networks in support of goal maintenance. Nat
Neurosci. 2015 Sep;18(9):1318-24. PubMed Central PMCID: PMC4551604.
d. Voytek B, D'Esposito M, Crone N, Knight RT. A method for event-related phase/amplitude
coupling. Neuroimage. 2013 Jan 1;64:416-24. PubMed Central PMCID: PMC3508071.

The neural basis of age-related cognitive decline: Drawing on my past work showing that aperiodic
“noise” is associated with cognitive impairments and age-related working memory decline, my lab has
begun work exploring the more precise neural origins of aperiodic activity and how it interacts with
oscillations and cognition. We tested the role that oscillatory disruptions contribute to age-related
cognitive decline by showing that oscillations in visual cortex in healthy older adults fail to reset in
response to informative cues that alert them to upcoming working memory challenge (a). Importantly,
this failure of oscillatory disruptions predicted the degree of working memory impairment among older
adults. In another study, we show that the increased putative EI “noise” among older adults is linked to
oscillatory disruptions, which in turn are linked to visual attention impairments among older adults (b).
Using our novel analytical framework for separating oscillatory and aperiodic activity, funded by our
current NIGMS R01, we showed that both oscillatory visual cortical alpha power and aperiodic power
predict individual visual working memory performance among older adults (c). Finally, also with our
current NIGMS R01 funding, we have also adopted a "data science" approach of combining multiple large, heterogenous human brain datasets to show that this aperiodic neural activity is mathematically
related to population neuronal timescales. Whole-brain iEEG mapping of these timescales shows a
gradient that increases along the principal sensorimotor-to-association axis across the entire human
cortex, which shows direct alignment with the expression of excitation- and inhibition-related genes,
as well as genes specific to voltage-gated transmembrane ion transporters. Finally, we showed that
neuronal timescales are functionally dynamic: prefrontal cortex timescales expand during working
memory maintenance and predict individual performance, while cortex-wide timescales compress with
aging (d).
a. Donoghue T, Haller M, Peterson EJ, Varma P, Sebastian P, Gao R, Noto T, Lara AH, Wallis JD,
Knight RT, Shestyuk A, Voytek B. Parameterizing neural power spectra into periodic and aperiodic
components. Nat Neurosci. 2020 Dec;23(12):1655-1665. PubMed Central PMCID: PMC8106550.
b. Gao R, van den Brink RL, Pfeffer T, Voytek B. Neuronal timescales are functionally dynamic and
shaped by cortical microarchitecture. Elife. 2020 Nov 23;9 PubMed Central PMCID: PMC7755395.
c. Tran TT, Rolle CE, Gazzaley A, Voytek B. Linked Sources of Neural Noise Contribute to Age-
related Cognitive Decline. J Cogn Neurosci. 2020 Sep;32(9):1813-1822. PubMed Central PMCID:
PMC7474516.
d. Tran TT, Hoffner NC, LaHue SC, Tseng L, Voytek B. Alpha phase dynamics predict age-related
visual working memory decline. Neuroimage. 2016 Dec;143:196-203. PubMed PMID: 27577719.


## My papers (titles + abstracts)
* Voytek JB & Voytek B. Automated cognome construction and semi-automated hypothesis generation. J Neurosci Methods 208, 92-100 (2012).

* Voytek B, Kayser A, Badre D, Fegen D, Chang EF, Crone NE, Parvizi J, Knight RT, D’Esposito M. Oscillatory dynamics coordinating human frontal networks in support of goal maintenance. Nature Neurosci 18, 1318-1324 (2015).

* Voytek B, Kramer MA, Case J, Lepage KQ, Tempesta ZR, Knight RT, Gazzaley A. Age-related changes in 1/f neural electrophysiological noise. J Neurosci 35, 13257-13265 (2015).

* Tran TT, Hoffner NC, LaHue SC, Tseng L, Voytek B. Alpha phase dynamics predict age-related visual working memory decline. NeuroImage 143, 196-203 (2016).

* Cole SR, van der Meij R, Peterson EJ, de Hemptinne C, Starr PA, Voytek B. Nonsinusoidal beta oscillations reflect cortical pathophysiology in Parkinson's disease. J Neurosci 37, 4830-4840 (2017).

* Voytek B, Samaha J, Rolle CE, Greenberg Z, Gill N, Porat S, Kader T, Rahman S, Malzyner R, Gazzaley A. Preparatory encoding of the fine scale of human spatial attention. J Cogn Neurosci 29, 1302-1310 (2017).

* Gao RD, Peterson EJ, Voytek B. Inferring synaptic excitation/inhibition balance from field potentials. NeuroImage 158, 70-78 (2017).

* van der Meij R & Voytek B. Uncovering neuronal networks defined by consistent between-neuron spike timing from neuronal spike recordings. eNeuro 5, (2018).

* Moore SM, Seidman JS, Ellegood J, Gao R, Savchenko A, Troutman TD, Abe Y, Stender J, Lee D, Wang S, Voytek B, Lerch JP, Suh H, Glass C, Muotri A. Setd5 haploinsufficiency alters neuronal network connectivity and leads to autistic-like behaviors in mice. Translational Psychiatry 9, 1-12 (2019).

* Jackson N, Cole SR, Voytek B, Swann NC. Characteristics of waveform shape in Parkinson’s disease detected with scalp electroencephalography. eNeuro 6, (2019).

* Veerakumar A, Tiruvadi V, Howell B, Waters AC, Crowell AL, Voytek B, Posse PR, Denison L, Rajendra JK, Edwards JA, Bijanki KR, Choi KS, Mayberg HS. Field potential 1/f activity in the subcallosal cingulate region as a candidate signal for monitoring deep brain stimulation for treatment resistant depression. J Neurophysiol 122, 1023-1025 (2019).

* Cole S & Voytek B. Cycle-by-cycle analysis of neural oscillations. J Neurophysiol 122, 849-861 (2019).

* Trujillo CA*, Gao R*, Negraes PD*, Chaim IA, Domissy A, Vandenberghe M, Devor A, Yeo GW, Voytek B#, Muotri AR#. Complex Oscillatory Waves Emerging from Cortical Organoids Model Early Human Brain Network Development. Cell Stem Cell 25, 558-569 (2019). *,#equal contribution

* Robertson MM, Furlong S, Voytek B, Donoghue T, Boettiger CA, Sheridan MA. EEG Power Spectral Slope differs by ADHD status and stimulant medication exposure in early childhood. J Neurophysiol 122, 2427-2437 (2019).

* Molina JL, Voytek B, Thomas ML, Joshi YB, Bhakta SG, Talledo JA, Swerdlow NR, Light GA. Memantine effects on EEG measures of putative excitatory/inhibitory balance in schizophrenia. Biol Psychiatry Cogn Neurosci Neuroimaging 5(6), 562-568 (2020).

* Ghatak S, Dolatabadi N, Gao R, Wu Y, Scott H, Trudler D, Sultan A, Ambasudhan R, Nakamura T, Masliah E, Talantova M, Voytek B, Lipton SA. NitroSynapsin ameliorates hypersynchronous neural network activity in Alzheimer hiPSC models. Mol Psychiatry (2020).

* Tran TT, Rolle CE, Gazzaley A, Voytek B. Linked sources of neural noise contribute to age-related cognitive decline. J Cogn Neurosci 32(9), 1813-1822 (2020).

* Peterson EJ & Voytek B. Homeostatic mechanisms may shape the type and duration of oscillatory modulation. J Neurophysiol 124, 168-177 (2020).

* Donoghue T, Dominguez J, Voytek B. Electrophysiological frequency band ratio measures conflate periodic and aperiodic neural activity. eNeuro (2020).

* Donoghue T*, Haller M*, Peterson E*, Varma P, Sebastian P, Gao R, Noto T, Lara AH, Wallis JD, Knight RT, Shestyuk A#, Voytek B#. Parameterizing neural power spectra into periodic and aperiodic components. Nature Neurosci 23, 1655-65 (2020). *,#equal contribution

* Gao R, van den Brink RL, Pfeffer T, Voytek B. Neuronal timescales are functionally dynamic and shaped by cortical microarchitecture. eLife e61277 (2020).

* Schaworonkow N & Voytek B. Longitudinal changes in aperiodic and periodic activity in electrophysiological recordings in the first seven months of life. Dev Cogn Neurosci (2021).

* Donoghue T, Voytek B, Ellis S. Teaching creative and practical data science at scale. J Stat Data Sci Educ 29, S27-39 (2021).

* Brown DE, Chavez JI, Nguyen DH, Kadwory A, Voytek B, Arneodo E, Gentner TQ, Gilja V. Local field potentials in a pre-motor region predict learned vocal sequences. PLOS Comput Biol (2021).

* Schaworonkow N & Voytek B. Enhancing oscillations in intracranial electrophysiological recordings with data-driven spatial filters. PLOS Comput Biol (2021).

* Waschke L, Donoghue T, Fiedler L, Smith S, Garrett DD, Voytek B#, Obleser J#. Modality-specific tracking of attention and sensory statistics in the human electrophysiological spectral exponent. eLife e70068 (2021). #equal contribution

* Donoghue T & Voytek B. Automated meta-analysis of the event-related potential (ERP) literature. Sci Rep (2022).

* Ostlund B, Donoghue T, Anaya B, Gunther KE, Karalunas SL, Voytek B, Pérez-Edgar KE. Spectral parameterization for studying neurodevelopment: How and why. Dev Cogn Neurosci (2022).

* Donoghue T, Voytek B, Ellis S. Course materials for data science in practice. J Open Source Educ (2022).
Peterson E, Rosen B, Belger A, Voytek B, Campbell A. Aperiodic neural activity is a better predictor of schizophrenia than neural oscillations. Clin EEG Neurosci (2023).

* Shafiei G, Fulcher BD, Voytek B, Satterthwaite TD, Bailler S, Misic B. Neurophysiological signatures of cortical micro-architecture. Nature Comms (2023).

* Smith SE, Ma V, Gonzalez C, Chapman A, Printz D, Voytek B#, Soltani M#. Clinical EEG slowing induced by electroconvulsive therapy is better described by increased frontal aperiodic activity. Transl Psychiatry (2023). #equal contribution

* Smith SE*, Kosik EL*, van Engen Q*, Hill AT, Zomorrodi R, Daskalakis ZJ, Hadas I#, Voytek B#. Magnetic seizure therapy and electroconvulsive therapy increase aperiodic activity. Transl Psychiatry (2023). *,#equal contribution

* Martin-Burgos B, McPherson T, Hammonds R, Gao R, Muotri A, Voytek B. Development of neuronal timescales in human cortical organoids and rat hippocampus dissociated cultures. J Neurophysiol (2024).

* Monchy N*, Modolo J*, Houvenaghel JF, Voytek B#, Duprez J#. Changes in electrophysiological aperiodic activity during cognitive control in Parkinson’s disease. Brain Communications (2024). *,#equal contribution

* Bender A, Voytek B#, Schaworonkow N#. Resting-state alpha and mu rhythms change shape across development but lack diagnostic sensitivity for ADHD and autism. J Cogn Neurosci (2025). #equal contribution

* Preston MJ, Schaworonkow N, Voytek B. Time-resolved aperiodic and oscillatory dynamics during human visual memory encoding. J Neurosci (2025).

* Bender A, Zhao C, Vogel E, Awh E, Voytek B. Differential representations of spatial location by aperiodic and alpha oscillatory activity in working memory. Proc Natl Acad Sci USA (2025).

* Monchy N*, Duprez J*, Houvenaghel JF, Legros A, Voytek B#, Modolo J#. Functional connectivity is dominated by aperiodic, rather than oscillatory, coupling. J Neurosci (2025). *,#equal contribution

* Cazares C, Hutton A, Paez G, Trauner D, Voytek B. Cannabidiol blood metabolite levels after cannabidiol treatment are associated with broadband EEG changes and improvements in visuomotor and non-verbal cognitive abilities in boys with autism requiring higher levels of support. Translational Psychiatry (2025).

* van Engen Q, Chau G, Smith A, Adam K, Donoghue T, Voytek B. Dissociating contributions of theta and alpha oscillations from aperiodic neural activity in human visual working memory. J Neurosci (2026).

* Pré D*, Cazares C*, WootenAT, Zhou H, Onofre I, Neil A, Logan T, Hu R, Lui JH, Voytek B#, Bang AG#. Pharmacological manipulation of nested oscillations in human iPSC-derived 2D neuronal networks. Neurobiol Aging(2026). *,#equal contribution

* Voytek B & Knight RT. Dynamic network communication as a unifying neural basis for cognition, development, aging, and disease. Biol Psychiatry 77, 1089-1097 (2015).

* Voytek B. The virtuous cycle of a data ecosystem. PLOS Comput Biol 12, 1-6 (2016).

* Cole SR & Voytek B. Brain oscillations and the importance of waveform shape. Trends Cogn Sci 21, 137-149 (2017).

* Voytek B. Social Media, Open Science, and Data Science are Inextricably Linked. Neuron 96, 1-4 (2017).

* Hanganu-Opatz I, Butt S, Hippenmeyer S, De Marco Garcia N, Cardin J, Voytek B, Muotri A. The logic of developing neocortical circuits in health and disease. J Neurosci (2021).

* Donoghue T, Schaworonkow N, Voytek B. Methodological considerations for studying neural oscillations. Eur J Neurosci (2021).

* Ahmad J, Ellis C, Leech R, Voytek B, Garces P, Jones E, Buitelaar J, Loth E, dos Santos F, Amil A, Verschure P, Murphy D, McAlonan G. From mechanisms to markers: novel non-invasive EEG proxy markers of the neural excitation and inhibition system in humans. Transl Psychiatry (2022).

* Voytek B. The data science future of neuroscience theory. Nature Methods (2022).

* Mushtaq F, …(84 authors)…, Valdes-Sosa P. 100 years of EEG for brain and behaviour research. Nature Hum Behav (2024)

* Fitzgerald M, Kosik E, Voytek B. Neuroscience: Unveiling hidden sources of noise. eLife (2024)
van Bree S, Levenstein D, Krause MR, Voytek B, Gao R. Processes and Measurements: a Framework for Understanding Neural Oscillations in Field Potentials. Trends Cogn Sci (2025)