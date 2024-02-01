# Publications

Below are publications from the sourmash team.

## sourmash fundamentals

[Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2), Irber et al., 2022. This is the core technical paper describing both FracMinHash and `sourmash gather`.

[Large-scale sequence comparisons with sourmash](https://f1000research.com/articles/8-1006),
Pierce et al., 2019. This is the original sourmash use case paper.

## Evaluation and benchmarking

[Evaluation of taxonomic classification and profiling methods for long-read shotgun metagenomic sequencing datasets](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05103-0),
Portik et al., 2022. This paper shows that sourmash is extremely
sensitive and very specific for taxonomic profiling.

## Petabase-scale search

[Biogeographic Distribution of Five Antarctic Cyanobacteria Using Large-Scale k-mer Searching with sourmash branchwater](https://www.biorxiv.org/content/10.1101/2022.10.27.514113v1),
Lumian et al., 2022. This paper uses sourmash and branchwater to
search ~500,000 public metagenomes for 5 query genomes, validates the
results using mapping, and discusses the biogeography of the query
species.

[Sourmash Branchwater Enables Lightweight Petabyte-Scale Sequence Search](https://www.biorxiv.org/content/10.1101/2022.11.02.514947v1),
Irber et al., 2022. This paper describes the technical underpinnings
of the first version of sourmash branchwater, for petabase-scale search.

## Advanced uses of sourmash

[Single-cell transcriptomics for the 99.9% of species without reference genomes](https://www.biorxiv.org/content/10.1101/2021.07.09.450799v1.abstract),
Botvinnik et al., 2021. This paper uses sourmash (and many other
techniques!) to analyze single cell data from the Chinese horseshoe
bat.

[Meta-analysis of metagenomes via machine learning and assembly graphs reveals strain switches in Crohn’s disease](https://www.biorxiv.org/content/10.1101/2022.06.30.498290v1.abstract),
Reiter et al., 2022. This paper uses sourmash and
[spacegraphcats](https://spacegraphcats.github.io/spacegraphcats/) to
detect and analyze strain-specific signals in fecal microbiomes from
the iHMP.

[Protein k-mers enable assembly-free microbial metapangenomics](https://www.biorxiv.org/content/10.1101/2022.06.27.497795v1),
Reiter et al., 2022. This paper develops a technique to use protein
k-mers to analyze metapangenome graphs from metagenomes.

## Additional works

Dr. Luiz Irber's PhD thesis,
[Decentralizing Indices for Genomic Data](https://github.com/luizirber/phd/releases),
describes several additional features of the sourmash ecosystem,
including wort, which monitors the SRA for new data sets and sketches
them automatically.
