# Using sourmash's built-in taxonomy handling from Python

sourmash works with taxonomies by connecting **identifiers** (typically
Genbank accessions or private identifiers) to **lineages** in lineage
databases.

Lineage databases can be loaded from CSV files or SQL databases; see
`sourmash tax prepare` to convert between the two.

For example,

~~~
>>> from sourmash.tax.tax_utils import MultiLineageDB, RankLineageInfo, LineagePair

>>> taxdb = MultiLineageDB.load(['doc/taxtax.csv'])
>>> print(f"Loaded {len(taxdb)} taxonomic identifiers.")
Loaded 3 taxonomic identifiers.

~~~

Accessions or identifiers are keys in the database:

~~~
>>> accession = 'GCF_014075335.1'
>>> accession in taxdb
True

~~~

Lineage database values are **tuples of `LineagePairs`**, or lineages -

~~~
>>> lineage = taxdb[accession]
>>> for linpair in lineage:
...    print(linpair.rank, linpair.name)
superkingdom d__Bacteria
phylum p__Proteobacteria
class c__Gammaproteobacteria
order o__Enterobacterales
family f__Enterobacteriaceae
genus g__Escherichia
species s__Escherichia coli

~~~

These lineages can be displayed, queried, and manipulated in a variety
of convenient ways by using `RankLineageInfo`:

~~~
>>> lineage = RankLineageInfo(lineage=lineage)
>>> lineage
RankLineageInfo(lineage_str='d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli')

~~~

Format the lineage for display with `display_lineage`:

~~~
>>> lineage.display_lineage()
'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli'

~~~

Confirm the presence (or absence) of specific ranks:

~~~
>>> lineage.rank_is_filled('species')
True
>>> lineage.rank_is_filled('strain')
False

~~~

Remove ranks below a particular rank with `pop_to_rank`:

~~~
>>> lin2 = lineage.pop_to_rank('class')
>>> lin2.display_lineage()
'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria'

~~~

If you want to get the name at a given rank, use `name_at_rank`:
~~~
>>> lin2.name_at_rank('phylum')
'p__Proteobacteria'

~~~

and you can get get the lineage tuple at a specific rank with
`lineage_at_rank`:

~~~
>>> lin2.lineage_at_rank('phylum')
(LineagePair(rank='superkingdom', name='d__Bacteria'), LineagePair(rank='phylum', name='p__Proteobacteria'))

~~~

## Calculating the lowest common ancestor of two or more lineages.

If you have multiple accessions, you can also find the lowest common
ancestor of their lineages; this functionality is used extensively in
both the `lca` and `tax` submodules of sourmash.

First, load the lineages:
~~~
>>> acclist = ['GCF_014075335.1', 'GCF_900013275.1', 'GCF_900698905.1']
>>> lineages = [ taxdb[acc] for acc in acclist ]
>>> lineages = [ RankLineageInfo(lineage=lin) for lin in lineages]
>>> for lin in lineages:
...    print(lin.display_lineage())
d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales_A;f__Enterobacteriaceae_A;g__Buchnera;s__Buchnera aphidicola_W

~~~

Then, calculate the LCA. Here, the LCA is just another `RankLineageInfo`
object that we update each time. Note that the LCA here is at class
`c__Gammaproteobacteria`, the lowest rank shared by all three lineages
(above `o__Enterobacterales` and `o__Enterobacterales_A`).

~~~
>>> lca = lineages.pop()
>>> while lineages:
...    lin = lineages.pop()
...    lca = lca.find_lca(lin)
>>> lca.display_lineage()
'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria'

~~~

Lineages can contain some or all of the ranks - you can get the
lowest rank in a lineage with `lowest_rank`,
~~~
>>> lca.lowest_rank
'class'

~~~

You can check if a rank is available in the taxonomy with
`check_rank_availability`, and check if that rank is filled with
`rank_is_filled`:

~~~
>>> lca.check_rank_availability('class')
True
>>> lca.check_rank_availability('family')
True
>>> lca.rank_is_filled('class')
True
>>> lca.rank_is_filled('family')
False

~~~

## Creating/initializing "traditional" lineages with `RankLineageInfo`

NCBI and GTDB taxonomies both use the seven ranks - superkingdom
(or domain), phylum, class, order, family, genus, and species - along
with an (optional) strain. Lineages from both taxonomies can be
loaded into and manipulated using `RankLineageInfo`.

You can load from a lineage string with semicolon separators,

~~~
>>> lin1 = RankLineageInfo(lineage_str='d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria')
>>> lin1
RankLineageInfo(lineage_str='d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria')

~~~

or build your own `LineagePair` tuples and supply them via `lineage`:
~~~
>>> lintup = (LineagePair(rank='superkingdom', name='d__Bacteria'),
...           LineagePair(rank='phylum', name='p__Proteobacteria'),
...           LineagePair(rank='class', name='c__Gammaproteobacteria'))
>>> lin2 = RankLineageInfo(lineage=lintup)
>>> lin2
RankLineageInfo(lineage_str='d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria')

~~~

Of course, lineages initialized either way are equivalent:

~~~
>>> lin1 == lin2
True

~~~

The taxonomy ranks themselves can be displayed with `taxlist` and
`ascending_taxlist` -

~~~
>>> lin1.taxlist
('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain')

>>> lin1.ascending_taxlist
('strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom')

~~~

~~~
>>> obj = RankLineageInfo(lineage_str='')
>>> type(obj)
<class 'sourmash.tax.tax_utils.RankLineageInfo'>
>>> obj
RankLineageInfo(lineage_str='')

~~~

## To discuss with Tessa:

1. Why does `taxdb[accession]` not return a `RankLineageInfo` object?
(And/or should it?)

2. what's the difference in use case for `lineage_at_rank` and `pop_to_rank`? How does this differ from pop_to_rank? It appears that it just gives you the tuple of LineagePairs rather than giving you a new RankLineageInfo object.

~~~
>>> lineage = RankLineageInfo(lineage_str='d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli')
>>> lineage.pop_to_rank('class')
RankLineageInfo(lineage_str='d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria')
>>> lineage.lineage_at_rank('class')
(LineagePair(rank='superkingdom', name='d__Bacteria'), LineagePair(rank='phylum', name='p__Proteobacteria'), LineagePair(rank='class', name='c__Gammaproteobacteria'))

~~~

`pop_to_rank` is used twice in the code base, `lineage_at_rank` is
used once.  Suggest picking one, or making a new method and using
that - maybe `lineage_to_rank`, which returns a RankLineageInfo
object?

3. suggest changing `check_rank_availability`

`check_rank_availability` is all about raising an exception, it seems -
no caller checks its return value. Is it OK to remove the return value?

## TODO items @CTB:

- [ ] discuss/describe NCBI taxids, zip_taxid, display_taxid
- [ ] describe zip_lineage
- [ ] document/test is_compatible, is_lineage_match

