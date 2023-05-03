# wrappers

Snakemake wrappers adhering to Snakelines focus on less technical users.

This repository is an alternative to the official [snakemake wrappers](https://github.com/snakemake/snakemake-wrappers). Official wrappers focus on more technical users, exposing in a rule usually only the `extra` parameter which summarizes all arguments and parameters. For our use case, wrappers should expose parameters individually thus removing the responsibility of the user to pass the correct parameter values to the wrapper and instead placing the responsibility to the developers' shoulders. Many wrappers are however simple having no parameters, in that case should be used from the official repository.

## Development

Install snakemake, for example in the environment `snakemake_dev`:

```shell
mamba create -c conda-forge -c bioconda --name snakemake_dev snakemake pre_commit
```

Then set up `pre-commit` in the repository:

```bash
pre-commit install
```

Now before any commit, a defined set of actions will be performed, such as linting, formatting, etc.

### Wrapper directory structure

Similar structure to the official snakemake wrappers is used. Each wrapper must be placed in the directory:

```bash
wrappers/{tool}/{command}
```

Where `tool` is the name of the wrapped tool and `command` is the name of the command, for example `wrappers/bwa/index`.

Each wrapper must adhere to the following structure:

```bash
README.md
environment.yaml
test/Snakefile
test/input
wrapper.py
```

- `README.md` - here provide a simple description of the wrapped tool.
- `environment.yaml` - the conda environment file used by the wrapper.
- `wrapper.py` - the wrapper logic. Here you have access to the `snakemake` object passed from rules calling the wrapper, e.g. access to input, output, etc. by `snakemake.input[0]` for example.
- test/Snakefile - an example of a rule calling the wrapper.
- test/input - test inputs for the example Snakefile.

Adhering to the structure is important, as github actions rely on this structure to automatically run tests, build documentation, etc.

### Commits and PRs

In commits you should follow the [Conventional Commits spec](https://www.conventionalcommits.org/en/v1.0.0/). Each PR is automatically validated by the GH action.

Further, any push (i.e. after merged PR) to the `main` branch results in a new PR with:

- a new release following the [Semantic Versioning](https://semver.org/)
- an automatic changelog as parsed from the commit history

## Wrapper example

Here we provide a simple example for the bwa index wrapper. It is stored in the `wrapper/bwa/index` directory.

Example of the structure at `wrapper/bwa/index`:

```bash
README.md
environment.yaml
test/Snakefile
test/input/ref.fa
wrapper.py
```

Example of the Snakefile:

```python
rule bwa__build_index:
    input:
        "ref.fa",
    output:
        idx=multiext("bwa_index/ref", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output.idx[0])[0],
        approach="bwtsw",
    log:
        "logs/bwa/build_index.log",
    wrapper:
        "https://github.com/xsitarcik/wrappers/raw/main/wrappers/bwa/index"
```

Example of the wrapper:

```python
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("bwa index -p {snakemake.params.prefix} -a {snakemake.params.approach} {snakemake.input[0]} {log}")
```

Example of the environment:
```yaml
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - bwa=0.7.17
```

