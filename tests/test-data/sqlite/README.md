# test files for SqliteIndex etc. functionality

`prot.sqlmf` is a SQL version of the manifest in `tests/test-data/prot/all.zip`.

`delmont-6.csv` is a fixed-up version of `tests/test-data/lca/delmont-6.csv` that works with `sourmash tax`.

`lca.sqldb` is an `LCA_SqliteDatabase` created with `TARA_ASE_MAG_00031` and `TARA_PSW_MAG_00136`, using the lineage in `delmont-6.csv`.

`test.taxonomy.db` is a SqliteLineage v1.0 lineage db created with `sourmash tax prepare` from `tests/test-data/tax/test.taxonomy.db`.

`index.sqldb` is a k=31 sqldb created from `tests/test-data/{47,63}.fa.sig`.

`shewanella-lineage.csv` is a hand-hacked file containing lineages for 47 and 63.

`lca2.sqldb` is an `LCA_SqliteDatabase` created from `tests/test-data/{47,63}.fa.sig` and `shewanella-lineage.csv`.
