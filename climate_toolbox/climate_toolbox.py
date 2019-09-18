"""
This file describes the process for computing weighted climate data
"""

from __future__ import absolute_import

import climate_toolbox.io
import logging

logger = logging.getLogger('uploader')

def transform_and_weighted_aggregate_climate_data(
        climate_data_loader=None,
        loader_kwargs=None,
        transform=None,
        extra_transform_kwargs=None,
        aggregator=None,
        extra_aggregation_kwargs=None,
        validate=None,
        extra_validation_kwargs=None,
        writers=None,
        extra_writer_kwargs=None,
        metadata=None,
        interactive=False,
        logger=None,
        **iteration_kwargs):
    '''

    '''

    if climate_data_loader is None:
        climate_data_loader = (
            climate_toolbox.io.load_and_standardize_climate_dataset_from_pattern)

    if loader_kwargs is None:
        loader_kwargs = {}

    if transform is None:
        transform = lambda x, **args, **kwargs: x

    if extra_transform_kwargs is None:
        extra_transform_kwargs = {}

    if aggregator is None:
        aggregator = lambda x, **args, **kwargs: x

    if extra_aggregation_kwargs is None:
        extra_aggregation_kwargs = {}

    if validator is None:
        validator = lambda x: **args, **kwargs, None

    if extra_validation_kwargs is None:
        extra_validation_kwargs = {}

    if writers is None:
        writers = [lambda  x, **args, **kwargs: None]

    if extra_writer_kwargs is None:
        extra_writer_kwargs = {}

    logger.debug('loading climate data')
    data = climate_data_loader(**loader_kwargs, **iteration_kwargs)

    logger.debug('transforming climate data')
    data = transform(data, **extra_transform_kwargs, **iteration_kwargs)

    logger.debug('aggregating transformed data')
    data = aggregator(data, **extra_aggregation_kwargs, **iteration_kwargs)

    logger.debug('validating results')
    data = validate(data, **extra_validation_kwargs, **iteration_kwargs)

    logger.debug('updating metadata')
    if metadata is not None:
        data.attrs.update(metadata, **iteration_kwargs)

    if interactive:
        return data

    logger.debug('writing results')
    for writer in writers:
        writer(data, **extra_writer_kwargs, **iteration_kwargs)

    logger.debug('done')
