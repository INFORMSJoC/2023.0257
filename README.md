![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)

# RMA: Ranking based on model averaging

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](https://github.com/INFORMSJoC/2023.0257/blob/master/LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [RMA: Ranking based on model averaging](https://doi.org/10.1287/ijoc.2023.0257) by Z. Feng, B. He, T. Xie, X. Zhang, and P. Xian. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0257

https://doi.org/10.1287/ijoc.2023.0257.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{FengRMA2024,
  author =        {Z. Feng, B. He, T. Xie, X. Zhang, and P. Xian},
  publisher =     {INFORMS Journal on Computing},
  title =         {RMA: Ranking based on model averaging},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0257.cd},
  url =           {https://github.com/INFORMSJoC/2023.0257},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0257},
}  
```

## Description

 Ranking problems are commonly encountered in practical applications, including order priority ranking, wine
 quality ranking, and piston slap noise performance ranking. The responses of these ranking applications are
 often considered as continuous responses and there is uncertainty on which scoring function is used to model
 the responses. In this paper, we address the scoring function uncertainty of continuous response ranking
 problems by proposing a Ranking Model Averaging (RMA) method. With a set of candidate models varied
 by scoring functions, RMA assigns weights for each model determined by a K-fold cross-validation criterion
 based on pairwise loss. 

This project contains three folders: `data`, `results`, `src`.

- `data`: This folder includes the data used in the paper.
- `results`: This folder contains the results of the experiments.
- `src`: This folder contains the source code and the code for experiment comparison.

## Data

The  datasets used for the numerical study are available in the `data` directory.

## Replicating

To reproduce each result in the paper, please run the file with the corresponding case number in the "src" directory. For example, to generate the numerical results of weight in Simulation 1 of Section 5.1 (Figure 2), please run MATLAB `src/Sim1/Simu_1/Figure1_Figure2.m`.

See the README.md file in each folder for a detailed description.

## Ongoing Development

This code is being developed on an on-going basis at the author's [GitHub site](https://github.com/xpzong/2023.0257).
