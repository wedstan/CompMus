# Make the World Dance

What makes a piece of music danceable? This project explores that question using audio features of AI-generated tracks, first across a corpus of 108 tracks and then through a structural analysis of two individual tracks.

**Live dashboard:** https://wedstan.github.io/CompMus/

Coursework for *Computational Musicology*, University of Amsterdam.

## Summary of findings

1. **Tempo is not the main driver.** Faster tracks are not systematically more danceable; the most danceable track in the corpus runs at 63 BPM, and 100 of the 108 tracks are faster.
2. **Emotional features matter more.** The most danceable tracks score higher on arousal and valence and lower on instrumentalness.
3. **Structure matters.** Self-similarity matrices suggest that danceable tracks show clearer repetition and rhythmic regularity.
4. **AI vs non-AI.** A random forest separates AI-generated from non-AI tracks with about 65% cross-validated accuracy (vs. a 54% majority baseline), relying most on instrumentalness and arousal and least on tempo.

## Methods

- **Corpus level:** standardised features (danceability, arousal, valence, tempo, instrumentalness), clustered heatmaps of the top and bottom five tracks, and scatter plots with linear fits.
- **Classification:** random forest (`ranger`, 500 trees) evaluated with 5-fold cross-validation repeated 20 times; impurity-based feature importance.
- **Track level:** chromagrams, cepstrograms (MFCCs), chroma- and timbre-based self-similarity matrices, key-template chordograms and tempograms.

Built in R with `flexdashboard`, `tidyverse`, `recipes`, `ranger`, `heatmaply` and `plotly`.

## Repository structure

| File | Description |
|---|---|
| `index.Rmd` | Source of the dashboard |
| `index.html` | Rendered dashboard (served via GitHub Pages) |
| `styles.css` | Dashboard styling |
| `compmus.R` | Helper functions for chroma, MFCC, self-similarity and tempogram analysis (provided by the course) |
| `compmus2025.csv` | Corpus-level features for the class corpus |

## Reproducing

Open `CompMus.Rproj` in RStudio and knit `index.Rmd`. The track-level plots require the course's feature files in a local `features/` folder, which are not included in this repository.
