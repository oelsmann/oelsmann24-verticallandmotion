"""Execute the two versioned notebooks with this exact Python environment."""

import argparse
from pathlib import Path
import sys

from jupyter_client.kernelspec import KernelSpecManager
import nbformat
from nbclient import NotebookClient


class LocalEnvironment(KernelSpecManager):
    def get_kernel_spec(self, kernel_name):
        spec = super().get_kernel_spec(kernel_name)
        spec.argv = [
            sys.executable,
            "-m",
            "ipykernel_launcher",
            "-f",
            "{connection_file}",
        ]
        return spec


def main():
    names = ["01_input_data.ipynb", "02_projection_comparison.ipynb"]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--notebook", action="append", choices=names, help="Execute only this notebook"
    )
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    for name in args.notebook or names:
        path = directory / name
        nb = nbformat.read(path, as_version=4)
        client = NotebookClient(
            nb,
            timeout=900,
            kernel_name="python3",
            resources={"metadata": {"path": str(directory)}},
        )
        client.create_kernel_manager()
        client.km.kernel_spec_manager = LocalEnvironment()
        client.execute()
        nbformat.write(nb, path)
        print(f"Executed diagnostics/{name}", flush=True)


if __name__ == "__main__":
    main()
