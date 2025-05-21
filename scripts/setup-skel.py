#!/usr/bin/env python3

import argparse
import pathlib
import re
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


def setup_module(
    source_dir: pathlib.Path,
    dest_dir: pathlib.Path,
    as_reference: bool = False,
) -> None:
    """
    Copy a module from the source directory to the destination directory, preparing the module for use.
    If as_reference is True, the module is copied as a reference module.
    """

    shutil.copytree(source_dir, dest_dir, symlinks=True)
    # remove directories and files that are not needed
    shutil.rmtree(dest_dir / "setup")
    (dest_dir / ".gitignore").unlink(missing_ok=True)
    # remove all files with the .nb.html extension
    for file in dest_dir.glob("*.nb.html"):
        file.unlink()
    if as_reference:
        # remove all files with the -live.Rmd extension
        for file in dest_dir.glob("*-live.Rmd"):
            file.unlink()
    else:
        # remove Rmd files *if there is a matching -live version*
        for rmd in dest_dir.glob("*.Rmd"):
            if (dest_dir / (rmd.stem + "-live.Rmd")).is_file():
                rmd.unlink()


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
        default=pathlib.Path("/etc/skel"),
        help="Directory to be used as the skel directory. (default: '%(default)s')",
    )
    parser.add_argument(
        "-m",
        "--modules",
        type=str,
        help="Modules to be included in the training, comma separated: e.g., 'module1,module2'",
    )
    parser.add_argument(
        "--reference-modules",
        type=str,
        help="Modules to be included in the training as reference modules (with completed notebooks)",
    )
    args = parser.parse_args()

    if not args.base_dir.is_dir():
        exit(f"Base directory {args.base_dir} does not exist or is not a directory")
    if not args.modules:
        exit("At least one module is required")

    modules = re.split(r"[,;\s]+", args.modules)
    reference_modules = (
        re.split(r"[,;\s]+", args.reference_modules) if args.reference_modules else []
    )

    if any(m not in ALL_MODULES for m in modules):
        exit(f"Invalid module(s) specified. Available modules: {ALL_MODULES}")
    if any(m not in ALL_MODULES for m in reference_modules):
        exit(f"Invalid reference module(s) specified. Available modules: {ALL_MODULES}")

    # check that modules are in the base directory
    if any(not (args.base_dir / m).is_dir() for m in modules):
        exit(f"One or more modules do not exist in the base directory {args.base_dir}")

    # run the link-data.sh script in the base directory
    link_script = args.base_dir / "scripts/link-data.sh"
    if not link_script.is_file():
        exit(f"The required script {link_script} does not exist or is not a file")
    subprocess.run(["bash", link_script], check=True)

    # create the target directory
    args.skel_dir.mkdir(parents=True, exist_ok=True)
    # create the shared-data link
    (args.skel_dir / "shared-data").symlink_to("/shared/data")

    # make the base directory in the skel
    target_base = args.skel_dir / args.base_dir.name
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
    if args.reference_modules:
        for module in reference_modules:
            setup_module(
                args.base_dir / module,
                target_base / module,
                as_reference=True,
            )

    # remove miscellaneous files and directories
    for item in REMOVE_MISC:
        item_path = target_base / item
        if item_path.is_file():
            item_path.unlink()
        elif item_path.is_dir():
            shutil.rmtree(item_path)


if __name__ == "__main__":
    main()
