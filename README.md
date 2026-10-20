# Systems Molecular Biology – computer practicals (UCO, 2026-27)

Everything you need for the practicals runs in your web browser: R, RStudio, all
the packages, the command-line tools and the data. Nothing to install.

## Start (first session, about 5 minutes)

1. Create a free account at [github.com](https://github.com) if you do not have one.
2. On this page, click the green **Code** button → **Codespaces** → **Create codespace on main**.
3. Wait while the codespace is created (2–4 minutes the first time). The data and
   scripts are unpacked into the folder `~/BMSi`.
4. RStudio opens in a new browser tab. If it does not, open the **Ports** tab at the
   bottom of the window and click the globe icon next to **RStudio (8787)**.

## Every following session

Go to [github.com/codespaces](https://github.com/codespaces) and open **the same
codespace**. Do not create a new one: your results are kept in the one you already have.

## Where things are

| What | Where |
| --- | --- |
| Data of each practical | `~/BMSi/PC<n>_...` |
| Scripts | `~/BMSi/scripts` |
| Results | the `results` folder inside each practical's folder |
| Terminal (practicals 1 and 2) | RStudio → **Terminal** tab; `conda activate bmsi_pc1` (or `bmsi_pc2`) |

## Free usage: what to keep in mind

- A GitHub account includes 120 core-hours a month; this codespace has 2 cores, so
  that is **60 hours of use a month**. The practicals need about 12 hours in class plus
  your own work.
- A codespace stops by itself after 30 minutes without activity; stop it yourself when
  you finish (github.com/codespaces → **…** → **Stop codespace**).
- Download the results you want to keep (RStudio → Files → More → Export) and
  **delete the codespace when the course ends**.

---

## For the instructor

- `.devcontainer/build/`: the image (R 4.5, Bioconductor 3.22, packages, conda
  environments, RStudio Server). The workflow `.github/workflows/build-image.yml` builds
  it, publishes it to GHCR and test-runs every practical at codespace size (2 CPU, 8 GB);
  the results table is in the summary of each workflow run.
- `.devcontainer/devcontainer.json`: what students get (prebuilt image, 2-core machine).
- `data/`: the `ForStudents_*.zip` files, unpacked into `~/BMSi` when a codespace is
  created. Files above 100 MB are not accepted by GitHub.
- After the first build, make the package public once: GitHub profile → Packages →
  `bmsi-practicals` → Package settings → Change visibility → Public.
