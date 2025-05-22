#!/usr/bin/env python3

"""
This script sets up a skeleton directory for user setup for Data Lab training modules.

The main function is to copy the a set of modules (as specified in a JSON input file),
from the base directory to the destination directory, removing unnecessary files and directories.

When setting up the modules, this script will also create symlinks to shared data files
that are expected to be found in the /shared directory on the host system.

"""

import argparse
import json
import pathlib
import shutil
import subprocess

ALL_MODULES = {
    "intro-to-R-tidyverse",
    "RNA-seq",
    "scRNA-seq",
    "scRNA-seq-advanced",
    "pathway-analysis",
    "machine-learning",
}

BASE_FILES = [
    "CITATION.cff",
    "gitignore.user",
    "LICENSE.md",
    "module_structure_detail.png",
    "README.md",
]

REMOVE_MISC = [
    "scRNA-seq/directory_structure.txt",
    "intro-to-R-tidyverse/results",
    "scRNA-seq/data/glioblastoma",
]


def remove_items(directory: pathlib.Path, items: list) -> None:
    """
    Remove items (files and directories) from a directory.
    """
    for item in items:
        item_path = directory / item
        if item_path.is_file():
            item_path.unlink()
        elif item_path.is_dir():
            shutil.rmtree(item_path)


def setup_module(
    source_dir: pathlib.Path,
    dest_dir: pathlib.Path,
    as_reference: bool = False,
) -> None:
    """
    Copy a module from the source directory to the destination directory,
    preparing the module for use in training.
    Setup files are removed, as are rendered notebooks and completed Rmd files.

    If `as_reference=True`, the module the `-live.Rmd` files are removed,
    and the completed Rmd files are kept.
    """

    # copy the module, leaving symlinks intact
    shutil.copytree(source_dir, dest_dir, symlinks=True)
    # remove directories and files that are not needed for training
    remove_items(
        directory=dest_dir,
        items=[
            "setup",
            ".Rproj.user",
            ".Rhistory",
            ".gitignore",
        ],
    )

    # remove all files with the .nb.html extension
    for file in dest_dir.glob("*.nb.html"):
        file.unlink()
    if as_reference:
        # remove all files with the -live.Rmd extension
        for file in dest_dir.glob("*-live.Rmd"):
            file.unlink()
    else:
        # remove Rmd files *if there is a matching -live version*
        for file in dest_dir.glob("*.Rmd"):
            if (dest_dir / (file.stem + "-live.Rmd")).is_file():
                file.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a skel directory for user setup for CCDL projects"
    )
    parser.add_argument(
        "-b",
        "--base-dir",
        type=pathlib.Path,
        default=pathlib.Path("/etc/skel-templates/training-modules"),
        help="Directory containing the modules to be included in the training",
    )
    parser.add_argument(
        "-s",
        "--skel-dir",
        type=pathlib.Path,
        help="Directory to be used as the output skel directory. This will be created if it does not exist",
    )
    parser.add_argument(
        "-m",
        "--module-file",
        type=pathlib.Path,
        help="A JSON file containing the list of modules to be included in the training",
    )
    args = parser.parse_args()

    if not args.base_dir.is_dir():
        exit(f"Base directory {args.base_dir} does not exist or is not a directory")
    if not args.module_file.is_file() and not args.module_file.suffix == ".json":
        exit("A JSON module file must be provided.")

    module_data = json.loads(args.module_file.read_text())
    modules = module_data.get("modules", [])
    reference_modules = module_data.get("reference-modules", [])

    if any(m not in ALL_MODULES for m in modules):
        exit(f"Invalid module(s) specified. Available modules: {ALL_MODULES}")
    if any(m not in ALL_MODULES for m in reference_modules):
        exit(f"Invalid reference module(s) specified. Available modules: {ALL_MODULES}")

    # check that modules are in the base directory
    if any(not (args.base_dir / m).is_dir() for m in modules):
        exit(f"One or more modules do not exist in the base directory {args.base_dir}")

    # run the link-data.sh script in the base directory to create required symlinks
    link_script = args.base_dir / "scripts/link-data.sh"
    if not link_script.is_file():
        exit(f"The required script {link_script} does not exist or is not a file")
    subprocess.run(["bash", link_script], check=True)

    # create the target directory
    args.skel_dir.mkdir(parents=True, exist_ok=True)
    # create the shared-data link
    (args.skel_dir / "shared-data").symlink_to("/shared/data")

    # make the base directory in the skel
    target_base = args.skel_dir / args.base_dir.resolve().name
    target_base.mkdir(parents=True, exist_ok=False)
    # copy the user base files
    for file in BASE_FILES:
        shutil.copy(args.base_dir / file, target_base / file)
    # put the gitignore file in the correct place
    (target_base / "gitignore.user").rename(target_base / ".gitignore")

    # set up the modules
    for module in modules:
        setup_module(args.base_dir / module, target_base / module)
    # set up the reference modules

    for module in reference_modules:
        setup_module(args.base_dir / module, target_base / module, as_reference=True)

    # remove miscellaneous files and directories
    remove_items(target_base, REMOVE_MISC)


if __name__ == "__main__":
    main()
