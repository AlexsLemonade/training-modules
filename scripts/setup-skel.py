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

USER_BASE_FILES = [
    "CITATION.cff",
    "LICENSE.md",
    "modules_structure_detail.png",
    "README.md",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a skel directory for user setup for CCDL projects"
    )
    parser.add_argument(
        ["b", "base_dir"],
        type=pathlib.Path,
        default=pathlib.Path("/etc/skel-templates/training-modules"),
        help="Directory containing the modules to be included in the training",
    )
    parser.add_argument(
        ["s", "skel_dir"],
        type=pathlib.Path,
        default=pathlib.Path("/etc/skel"),
        help="Directory to be used as the skel directory. (default: '%(default)s')",
    )
    parser.add_argument(
        ["m", "modules"],
        type=str,
        help="Modules to be included in the training, comma separated: e.g., 'module1,module2'",
    )
    args = parser.parse_args()

    module_set = set(re.split(r"[,;\s]+", args.modules))
    if not args.base_dir.is_dir():
        exit(f"Base directory {args.base_dir} does not exist or is not a directory")
    if not args.training_tag:
        exit("Training tag is required")
    if not args.modules:
        exit("At least one module is required")
    if any(m not in ALL_MODULES for m in module_set):
        exit(f"Invalid module(s) specified. Available modules: {ALL_MODULES}")

    # check that modules are in the base directory
    if any(not (args.base_dir / m).is_dir() for m in module_set):
        exit(f"One or more modules do not exist in the base directory {args.base_dir}")

    # run the link-data.sh script in the base directory
    link_script = args.base_dir / "scripts/link-data.sh"
    if not link_script.is_file():
        exit(f"The required script {link_script} does not exist or is not a file")
    subprocess.run(["bash", link_script], check=True)

    # create the target directory
    args.skel_dir.mkdir(parents=True, exist_ok=False)
    # create the shared-data link
    (args.skel_dir / "shared-data").symlink_to("/shared/data")

    # make the base directory in the skel
    target_base = args.skel_dir / args.base_dir.name
    target_base.mkdir(parents=True, exist_ok=False)
    # copy the user base files
    for file in USER_BASE_FILES:
        shutil.copy(args.base_dir / file, target_base / file)
    # put the gitignore file in the correct place
    (target_base / "gitignore.user").rename(target_base / ".gitignore")

    # copy the modules
    for module in module_set:
        module_dir = target_base / module
        shutil.copytree(args.base_dir / module, module_dir, symlinks=True)
        # remove directories and files that are not needed
        shutil.rmtree(module_dir / "setup")
        (module_dir / ".gitignore").unlink(missing_ok=True)
        # remove all files with the .nb.html extension
        for file in module_dir.glob("*.nb.html"):
            file.unlink()
        # remove Rmd files *if there is a matching -live version*
        for rmd in module_dir.glob("*.Rmd"):
            if (module_dir / rmd.stem).with_suffix("-live.Rmd").is_file():
                rmd.unlink()


if __name__ == "__main__":
    main()
