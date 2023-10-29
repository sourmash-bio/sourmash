# Installing sourmash

This tutorial should run without modification on Linux or Mac OS X,
under [Miniforge](https://github.com/conda-forge/miniforge).

You'll need about 5 GB of free disk space.

## Install miniforge

If you don't have the `mamba` command installed, you'll need to install
[miniforge](https://github.com/conda-forge/miniforge#install).

On Linux, this should work:
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b

~/miniforge3/bin/mamba init

echo 'source ~/.bashrc' > ~/.bash_profile
source ~/.bash_profile
```
otherwise, follow the instructions [here](https://github.com/conda-forge/miniforge#install).

## Install sourmash

To install sourmash, create a new environment named `smash` and install sourmash:

```
mamba create -y -n smash sourmash
```

and then activate:
```
conda activate smash
```

You should now be able to use the `sourmash` command:

```
sourmash info
```

Voila!
