# Repository publication checklist

Before publishing the contents of `/data/scdrugflow`:

- [ ] Confirm `python -m pytest -q` passes after the final detection-fraction fix.
- [ ] Optionally rerun the RA Nextflow workflow and inspect `quality_summary.json` plus the bounded `top20_gene_quality.tsv` values.
- [ ] Review `git status --short`, then `git diff --cached --name-only` after staging.
- [ ] Confirm `.h5ad` files, `results/`, `work/`, `main.nf.before_*`, ZIP patches, `.nextflow*`, credentials and unrelated code are **not** staged.
- [ ] Confirm the name and open-source MIT licence in `LICENSE` are appropriate for the code you are releasing.
- [ ] Open the GitHub README to confirm the image and commands render correctly.
- [ ] Confirm GitHub Actions Python tests pass after the initial push.
- [ ] Review and be able to explain the pipeline steps, code, limitations and AI-assisted work before using it in a job application.

The documented run was performed by the repository user on a public dataset; no confidential thesis or institutional scripts should be published.
