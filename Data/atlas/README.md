# Atlas Provenance and Adaptations

## Source

The base CocoYeo143 parcellation and its associated region labels were obtained from [NeMo 2.1: Network Modification Tool](https://github.com/kjamison/nemo#readme). The upstream NeMo documentation describes CocoYeo143 as a 143-region parcellation combining 100 Schaefer cortical parcels, 16 FreeSurfer aseg subcortical regions, and 27 SUIT cerebellar regions.

## Project adaptations

The atlas resources were adapted for the spatial grids and laterality-normalized analyses used in this project:

1. Files containing `_grid` were resampled to the corresponding analysis grid with nearest-neighbor interpolation to preserve discrete ROI labels.
2. Files containing `_flipLR` are left-right flipped variants generated from the corresponding unflipped atlas.
3. The lookup table retains the ROI identities used to align atlas labels, regional features, and white matter disconnection matrices.

These adapted files are project derivatives and are not official NeMo distributions.

## Attribution and licensing

NeMo is distributed under the MIT License and is copyright Keith W. Jamison. Users should consult the [upstream repository](https://github.com/kjamison/nemo) for the original software license, atlas documentation, and citations for the component parcellations.
