# Developer Notes

## Installation

### 1. Clone repository

```bash
git clone git@github.com:mlbendall/telescope.git telescope_dev.git
```

### 2. Create conda environment

An environment file is included in this repo. 

```bash
mamba env create -n telescope_dev -f telescope_dev.git/environment.yml
conda activate telescope_dev
```

### Checkout version (optional)

If not using the main branch, check out the correct version

```bash
cd telescope_dev.git
git pull #just in case anything changed
git checkout main
```

### Install using pip

Install `telescope` in interactive mode (`-e`) so that changes to the
code are reflected immediately without reinstalling. 

```bash
# change to repo directory if not already there
# cd telescope_dev.git

pip install -e . 
```


## Testing

The following one-liner will run telescope against the provided test
data and check the final log-likelihood calculation.

```bash
eval $(telescope test) 2>&1 | grep 'Final log-likelihood' | grep -q '95252.56293' && echo "Test OK" || echo "Test FAIL"
```
