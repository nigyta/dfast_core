#!/bin/env python
# coding: utf-8

import os
import platform
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


def set_binaries_path(app_root):
    os_name = platform.system()
    assert os_name == "Darwin" or os_name == "Linux"
    logger.info("OS type is {0}.".format(os_name))
    bin_dir = os.path.join(app_root, "bin", os_name)
    logger.info("Adding {0} to PATH.".format(bin_dir))
    os.environ["PATH"] = bin_dir + ":" + os.environ["PATH"]


if __name__ == '__main__':
    pass