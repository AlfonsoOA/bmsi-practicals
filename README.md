# Systems Molecular Biology – computer practicals (UCO, 2026-27)

Everything you need for the practicals runs in your web browser: R, RStudio, all
the packages, the command-line tools and the data. Nothing to install.

## Start (first session, about 5 minutes)

1. Create a free account at [github.com](https://github.com) if you do not have one.
2. **Before creating the codespace:** go to [github.com/settings/codespaces](https://github.com/settings/codespaces)
   → **Default idle timeout** → write **240** minutes → **Save**. Working in RStudio does
   not count as activity for GitHub, so with the default (30 minutes) the codespace
   stops in the middle of the practical. The setting only applies to codespaces created
   after changing it.
3. On this page, click the green **Code** button → **Codespaces** → **Create codespace on main**.
4. Wait while the codespace is created (2–4 minutes the first time). The data and
   scripts are unpacked into the folder `~/BMSi`.
5. When it is ready, the **Terminal** at the bottom of the window shows the address of
   RStudio (`https://<your-codespace>-8787.app.github.dev`): Ctrl + click on it. You can
   also open the **Ports** tab, right-click the **Forwarded address** of RStudio (8787) and
   choose **Open in Browser**. Never use `localhost:8787`: it does not work in the browser.

## Every session: one codespace, deleted at the end

1. At the start of the session create a new codespace (step 3 above; about 10 minutes until
   the Terminal says "RStudio is ready"). **Have only one codespace at a time.**
2. At the end, keep your results: `bash ~/BMSi/scripts/pack_results.sh PC1` (PC2, PC3, PC4),
   then RStudio → Files → `BMSi` → tick the zip → More → Export.
3. **Delete the codespace**: [github.com/codespaces](https://github.com/codespaces) → **…** → **Delete**.

## Where things are

| What | Where |
| --- | --- |
| Data of each practical | `~/BMSi/PC<n>_...` |
| Scripts | `~/BMSi/scripts` |
| Results | the `results` folder inside each practical's folder |
| Terminal (practicals 1 and 2) | RStudio → **Terminal** tab; `conda activate bmsi_pc1` (or `bmsi_pc2`) |

## Free usage: what to keep in mind

- A GitHub account includes 120 core-hours a month (this codespace has 4 cores, so **30 hours
  of use a month**, enough for the practicals) and **15 GB-month of storage**.
- **Storage is what runs out:** a codespace takes 10–15 GB for as long as it exists, even when
  stopped. Two weeks kept, or two codespaces at once, use up the month and GitHub blocks
  Codespaces until the next month. That is why each codespace is deleted at the end of the session.
- A codespace stops by itself after the idle timeout (step 2: 240 minutes); activity in the
  RStudio tab does not count, only in the codespace tab. If it stops in the middle of a session,
  nothing is lost: open it again from github.com/codespaces.

## For the instructor

- `.devcontainer/build/`: the image (R 4.5, Bioconductor 3.22, packages, conda
  environments, RStudio Server). The workflow `.github/workflows/build-image.yml` builds
  it, publishes it to GHCR and test-runs every practical at codespace size (2 CPU, 8 GB);
  the results table is in the summary of each workflow run.
- `.devcontainer/devcontainer.json`: what students get (prebuilt image, 4-core / 16 GB machine).
- `data/`: the `ForStudents_*.zip` files, unpacked into `~/BMSi` when a codespace is
  created. Files above 100 MB are not accepted by GitHub.
- After the first build, make the package public once: GitHub profile → Packages →
  `bmsi-practicals` → Package settings → Change visibility → Public.

## Licence and citation

Teaching materials: CC BY 4.0 (LICENSE-CONTENT.md). Code: MIT (LICENSE). How to cite: CITATION.cff.
