# rFRC4MINFLUX

This repository is distributed as an accompanying part of publication: [Weisong Zhao et al. Quantitatively mapping local quality at super-resolution scale by rolling Fourier ring correlation, <!-- Nature Biotechnology -->, X, XXX-XXX (2022)](https://www.nature.com/nbt/).

This repository is provided to reproduce our rFRC mapping on the MINFLUX data, including `Raw localization json files` provided by [Roman Schmidt et al.](https://doi.org/10.1038/s41467-021-21652-z); `Python scripts` to generate MINFLUX image stacks; `Archived version of rFRC ImageJ plugin (1/7 threshold)` to map local qualities; and `rFRC mapping results` for 3 representative cases.

[Portal](https://github.com/WeisongZhao/rFRC4MINFLUX/raw/main/rFRC_0143_-0.2.5.jar) to the `json files`.
[Portal](https://github.com/WeisongZhao/rFRC4MINFLUX/raw/main/rFRC_0143_-0.2.5.jar) to the `Python scripts`.
[Portal](https://github.com/WeisongZhao/rFRC4MINFLUX/raw/main/rFRC_0143_-0.2.5.jar) to the `ImageJ plugin`.
[Portal](https://github.com/WeisongZhao/rFRC4MINFLUX/raw/main/rFRC_0143_-0.2.5.jar) to the `rFRC mapping results`.


### Parameters for rFRC mapping 

For `STORM - Tubulin`:
- Block size: 64;
- Background intensity: 45;
- Skip: 1;
- Pixel size: 10 nm.

For `MINFLUX - Nup96`:
- Block size: 64;
- Background intensity: 60;
- Skip: 1;
- Pixel size: 1 nm.

For `MINFLUX - Beta II spectrin`:
- Block size: 64;
- Background intensity: 60;
- Skip: 1;
- Pixel size: 2 nm.

Because these MINFLUX examples are generated with a two-step filter, the image contains limited structural content, inducing large distance between the image-pair. Plus to avoid conservative resolution estimation, we used the 1/7 threshold in these representative cases. 

## Open source [rFRC4MINFLUX](https://github.com/WeisongZhao/rFRC4MINFLUX)
- This software and corresponding methods can only be used for **non-commercial** use, and they are under Open Data Commons Open Database License v1.0.
- Feedback, questions, bug reports and patches are welcome and encouraged!