# Bioconda recipe

Submitting this is a pull request against
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes), adding these files
under `recipes/printe/`.

Before opening that PR, fill in the `sha256` for the release tarball. Tag the release
first, then:

```bash
curl -sL https://github.com/cwb14/PrinTE/archive/refs/tags/v1.0.0.tar.gz | sha256sum
```

The placeholder in `meta.yaml` is all zeros, so a stale value fails loudly rather than
building the wrong source.

To try the recipe locally before submitting:

```bash
mamba install -c conda-forge conda-build
conda build conda/
```
