# Sequence Alignment

In Bioinformatics, a sequence alignment consists of rearranging sequences of primary structures, such as DNA, RNA or proteins, in order to maximize the similarity score. These similarities can be caused by functional, structural or even evolutionary relationships between them.
Sequence alignments can be divided in two categories:
- Global Alignment: the alignment is “forced” over the full extension of the sequences; - Local Alignment: the alignment tries to find similarity regions within long sequences
that might be overall really different.
The purpose of doing this is to homology: the similarity due to descent from a common ancestor. Often, can infer homology from similarity. If the two sequences originate from individuals with a common ancestor, mismatches can be interpreted as mutations, and gaps can be interpreted as insertions/deletions in one or another sequence that occurred ever since both species diverged in time.

## Test Genes
The chosen gene was HERC2 variations for Homo Sapiens (human) and Mus Musculus (domestic rat) species. HERC2 gene belongs to the HERC gene family that encodes unusually large proteins. Genetic variations are associated with skin/hair/eye pigmentation variability.

## Test Proteins
The chosen protein family was insulin-like growth factor-binding proteins (IGFBP) for Homo Sapiens (human), IGFBP-4, and Rattus Norvegicus (rat), IGFBP-5, species. “The IGFBPs help to lengthen the half-life of circulating insulin-like growth factors (IGFs) in all tissues, including the prostate. Individual IGFBPs may act to enhance or attenuate IGF signaling depending on their physiological context (i.e. cell type). Even with these similarities, some characteristics are different: chromosomal location, heparin binding domains, RGD recognition site, preference for binding IGF-I or IGF-II, and glycosylation and phosphorylation differences. These structural differences can have a tremendous impact on how the IGFBPs interact with cellular basement membranes.”

## Sequence Extraction
To extract gene sequences using gene IDs, I use the website https://www.ebi.ac.uk/ena as base.
To extract protein sequences using protein accession codes, I use the website www.ebi.ac.uk/proteins/api as base.
The returned info comes in fasta format. However, I format it in order to obtain the sequence alone. This way it works both with the server and with my algorithms.

## Scoring Matrices
PAM and BLOSUM are used as scoring matrices for proteins.
PAM matrices are highly criticized due to its assumption that each aminoacid in a sequence is equally mutable, which isn’t true. Another problem is that in the extrapolation of the PAM-1 matrix into higher order PAM-n matrices erros inherent in the PAM-1 matrix data are highly magnified. But the biggest issue is the dataset used to create PAM which basically consisted of globin proteins alone.
BLOSUM was created in order to address the issue of variable amino acid mutation rates within sequences. They are designed to improve the accuracy of alignments between distantly related protein sequences. Multiple alignments of related sequences were made, evolutionary distances were taken into account, etc. This is why BLOSUM is a better alternative to PAM.

## Server Connection
Using tools from the website https://www.ebi.ac.uk/ we are able to redeem a sequence and perform sequence alignment on both gene and protein sequences, while having the option to adjust the algorithms’ parameters.
The algorithms used using this resource are:

## Algorithms used:
- Global alignment algorithm: Needleman-Wunsch
- Local alignment algorithm: Smith-Waterman
- My implementations take into account affine cost. A gap of length k is more probable than k gaps of length 1; a gap may be due to a single mutational event while separated gaps are probably due to distinct mutational events. A linear gap penalty function treats these cases the same, so, in order to implement affine function, other than the 3 matrices strategy, I have h as a penalty associated with opening a gap, and g as a smaller penalty for extending the gap.



