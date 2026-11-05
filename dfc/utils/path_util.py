#!/bin/env python
# coding: utf-8

import os
import shutil
from logging import getLogger

logger = getLogger(__name__)


def create_output_directory(config):
    directory = config.WORK_DIR
    force_overwrite = config.FORCE_OVERWRITE
    if os.path.exists(directory):
        if os.path.isdir(directory):
            if force_overwrite:
                shutil.rmtree(directory)
                logger.warning("Working directory is set to an existing directory: '{}'.".format(directory))
                os.makedirs(directory)
            else:
                logger.error("Output directory '{}' already exists. To overwrite the existing directory, use '--force' option. Aborting...".format(directory))
                exit(1)
        else:
            logger.error("Cannot create output directory '{}'. Aborting...".format(directory))
            exit(1)
    else:
        try:
            os.makedirs(directory)
            logger.debug("Creating an output directory '{}'.".format(directory))

        except OSError:
            logger.error("Failed to create an output directory '{}'. Aborting...".format(directory))
            exit(1)


if __name__ == '__main__':
    pass