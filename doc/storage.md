# Storing SBTs

We support different storage options for the internal SBT data.


## Available storages

### FSStorage

The initial storage schema.
Saves internal SBT data in a hidden directory near the SBT JSON description.
- Pros: easy to create
- Cons: annoying to distribute (thousands of files).
        We used to create a tar file of JSON + hidden directory,
        which requires extracting and using more disk space.

### ZipStorage

Similar to FSStorage,
but saves the internal SBT data in a `zip` file.
- Pros
  * easy to distribute (one file)
- Cons
  * still need to distribute and download everything
    (you need the full zip file available locally to be able to use the SBT).

### IPFSStorage

Uses IPFS to store internal SBT data,
allowing partial database download.
- Pros
  * easy to distribute (one file, the SBT JSON description)
  * only data needed for analysis is downloaded
  * benefits from more people storing the data in their computers and sharing bandwidth
- Cons:
  * needs IPFS daemon running in the computer
  * takes longer to run if data is not prefetched

### RedisStorage

Meant to be a fast in-memory storage.
There won't be a public Redis server to provide the internal SBT data,
but this storage is a good option for loading data from others sources and sharing
with other processes or servers in your private network.
- Pros
  * Shareable between processes or servers in a network
  * Faster access time than reading from disk (probably?)
- Cons
  * No public server for the data (need to convert from other sources)

## Converting an existing tree to use a new storage

You can convert SBTs to another storage using the `sourmash storage convert` command:
``` bash
$ sourmash storage convert -b new_storage_type database.sbt.json
```

For example: to convert a tree to IPFSStorage, do
``` bash
$ sourmash storage convert -b ipfs database.sbt.json
```
