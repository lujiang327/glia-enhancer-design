#!/usr/bin/env python3
"""Run upstream ChromBPNet hyperparameter preparation with a fixed RNG seed."""
import numpy as np

from chrombpnet.helpers.hyperparameters import find_chrombpnet_hyperparams


if __name__ == '__main__':
    parser = find_chrombpnet_hyperparams.parse_data_args()
    parser.add_argument(
        '--sampling-seed', type=int, default=42,
        help='Seed for the upstream nonpeak subsample used to estimate count thresholds.'
    )
    args = find_chrombpnet_hyperparams.parse_model_args(parser)
    sampling_seed = args.sampling_seed
    del args.sampling_seed
    np.random.seed(sampling_seed)
    print(f'ChromBPNet hyperparameter sampling seed: {sampling_seed}')
    find_chrombpnet_hyperparams.main(args)
