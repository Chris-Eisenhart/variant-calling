# variant_calling
## A repo to explore variant calling from fastq files as well as some supporting infrastructure
Variant calling is the process of determining which genomic positions differ in a given sample from the reference genome.
This is typically done through a two step process; aligning the fastq reads to the reference genome, then calling variants
from the alignment files.

The repo uses BWA to perform alignments, then python code to read in the alignments and call variants.

## Setup
```sh
# Install dependencies
pipenv install --dev

# Setup pre-commit and pre-push hooks
pipenv run pre-commit install -t pre-commit
pipenv run pre-commit install -t pre-push
```

## Contributing
The style is black + isort + flake8, additionally type hinting is enforced via mypy. 

Keep all functions under 100 lines, take time to name variables appropriately, do not solve a previously
solved problem a different way.

Reach out to Chris to get invited to the JIRA board! 


