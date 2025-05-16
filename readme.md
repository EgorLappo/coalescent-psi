Code accompanying the article "Coalescent theory of the ψ directionality index" by Egor Lappo and Noah A. Rosenberg.

We use [`pixi`](https://pixi.sh/latest/) to manage the programming environment.
To run the code in this repository, first install `pixi` and then use `pixi run` to execute scripts.

* The file `plots.ipynb` generates Figures 4, 5, and 6 in the paper. To run it you can use `pixi run jupyter notebook`.

* The scripts in the `dmel_data` folder perform the bootstrapped computation of ψ for *Drosophila melanogaster* data (Section 6 of the paper).
  Just execute all numbered scripts in order: `00_download_data.py` downloads data from EMBL, `01_read_accessions.py` parses EMBL files, etc.
  The last script `05_aggregate.py` prints out the final values of mean and variance of ψ.

