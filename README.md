# Jupyter Notebook

## Installation

### Poetry

```sh
curl -sSL https://raw.githubusercontent.com/sdispater/poetry/master/get-poetry.py | python
export PATH=$HOME/.poetry/bin:$PATH

poetry completions zsh > $(brew --prefix)/share/zsh/site-functions/_poetry
fpath+=~/.zfunc

poetry install
```

### JupyterLab

```sh
jupyter lab
```
