#! /usr/bin/env python
# coding: UTF8

from .base_tools import Tool


class DefenseFinder(Tool):
    """DefenseFinder (anti-phage defense and anti-defense systems)

    Tool type: functional annotation (used by components.DefenseFinder)
    URL: https://github.com/mdmparis/defense-finder
    REF: Tesson et al., 2022
    """
    version = None
    NAME = "DefenseFinder"
    VERSION_CHECK_CMD = ["defense-finder", "version", "2>/dev/null"]
    VERSION_PATTERN = r"DefenseFinder version (\d+\.\d+\.\d+)"
    VERSION_ERROR_MSG = "DefenseFinder is required for --modules. Install it with 'conda install -c bioconda defense-finder'."

    def get_command(self, query_file, out_dir, models_dir, cpu, cmd_options=""):
        cmd = ["defense-finder", "run", "--db-type", "gembase", "--models-dir", models_dir,
               "-w", str(cpu), "-o", out_dir]
        if cmd_options:
            cmd += cmd_options.split()
        return cmd + [query_file]
