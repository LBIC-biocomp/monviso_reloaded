import subprocess
from pathlib import Path
from typing import Union


class Cobalt:
    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def run(
        self,
        hits_file: Union[str, Path],
        aligned_file: Union[str, Path],
        cobalt_home: str,
    ) -> bool:
        """_summary_

        Args:
            hits_file (Union[str,Path]): Path to the blastp output
            aligned_file (Union[str,Path]): Path where to direct the
            Cobalt output.

            cobalt_home (str): cobalt home, containing the executables

        Returns:
            bool: True if successful, else False
        """
        hits_file = Path(hits_file)
        aligned_file = Path(aligned_file)

        print("Doing MSA with COBALT.")
        command = [
            str(Path(cobalt_home,"cobalt")),
            "-i",
            str(hits_file),
            "-outfmt",
            "mfasta",
            "-end_gapopen",
            "5",
            "-end_gapextend",
            "1",
            "-gapopen",
            "11",
            "-gapextend",
            "1",
            "-blast_evalue",
            "0.003",
            "-norps",
            "T",
            "-treemethod",
            "clust",
        ]

        with aligned_file.open("w") as output_file:
            result = subprocess.run(
                command, stdout=output_file, stderr=subprocess.PIPE, text=True
            )

        if result.returncode == 0:
            print("MSA with COBALT completed successfully.")
            return True
        else:
            print(f"COBALT command failed with error: {result.stderr}")
            return False
