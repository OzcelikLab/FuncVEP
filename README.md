# FuncVEP

FuncVEP is a family of LightGBM classifiers for predicting the functional effect (Damaging / Neutral) of missense variants. It is trained on balanced and diverse functional data, enabling improved accuracy and generalization.

- FuncVEP Web Portal: https://funcvep.bilkent.edu.tr
- Precomputed FuncVEP scores for all possible missense variants (preprint version): https://zenodo.org/records/17036008
  - Updated scores will be made available upon publication.

> **Note:** This is the updated implementation of the FuncVEP framework. The code corresponding to the preprint is available in the `preprint` branch for reproducibility.

![FuncVEP Framework](resources/framework.png)

---

## Getting Started

### Requirements

You can install required Python packages using one of the two methods below:

#### Option 1: pip

```bash
pip install -r requirements.txt
```

#### Option 2: Conda

```bash
conda env create -f environment.yml
conda activate funcvep
```

---

## Trained Models

The trained models used in this project are included in the `models/` directory.

---

## Reproducing the Results

The project code is organized into numbered Jupyter notebooks that reproduce the full FuncVEP framework.

### Jupyter Notebooks (`notebooks/`)

Run the following notebooks in order:

1. `00_download_data.ipynb` (Downloads and extracts required data from Zenodo)
2. `01_feature_imputation.ipynb`
3. `02_train_test_set_construction.ipynb`
4. `03_model_training.ipynb`
5. `04_core_variant_scoring.ipynb`
6. `05_inference_on_new_variants.ipynb`
7. `06_model_benchmarking.ipynb`
8. `07_shap_analysis.ipynb`
9. `08_benchmark_figure_generation.ipynb`

---

## Citation

The FuncVEP paper will be provided as soon as it is available.

---

## License

This project is available under the PolyForm Strict License.

---
