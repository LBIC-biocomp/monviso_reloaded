import logging
import shutil
import os
import time
from pathlib import Path
import subprocess
from Bio import SeqIO
from Bio.Blast import NCBIWWW as blastq
from Bio.Blast import NCBIXML as blastparser


def run_cobalt(gene: str, cobalt_home: str) -> bool:
    """
    Runs the COBALT multiple sequence alignment tool on a specified gene.

    Parameters:
    - gene (str): The name of the gene to align.
    - cobalt_home (str): The path to the COBALT installation directory.

    Returns:
    - bool: True if COBALT ran successfully, False otherwise.
    """
    hits_file = Path(f"{gene}_hits.fasta")
    aligned_file = Path("aligned.fasta")

    if hits_file.exists() and hits_file.stat().st_size > 0:
        logging.info("Doing MSA with COBALT.")
        command = [
            f"{cobalt_home}/cobalt",
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
            logging.info("MSA with COBALT completed successfully.")
            return True
        else:
            logging.error(f"COBALT command failed with error: {result.stderr}")
            return False
    else:
        logging.warning(
            f"No hits file found or file is empty for gene: {gene}"
        )
        return False

def make_isof_folders(gene_path, iso: str) -> Path:
    """
    Creates a directory for a given isoform within a specified gene path.
    If the directory already exists, it will be emptied.

    Parameters:
    - gene_path (str or Path): The base directory path for the gene.
    - iso (str): The isoform identifier used to name the directory.

    Returns:
    - Path: The path to the isoform directory.
    """

    # Convert gene_path to Path object if it's not already one
    gene_path = Path(gene_path)
    outdir = gene_path / iso

    if outdir.exists():
        logging.warning(f"Directory {outdir} already exists. Emptying it.")
        for item in outdir.iterdir():
            if item.is_dir():
                # Recursively remove directories
                shutil.rmtree(item)
            else:
                # Remove files
                item.unlink()
    else:
        # Make the output directory if it doesn't exist
        outdir.mkdir(parents=True, exist_ok=True)

    return outdir


def start_query(myfasta, gene):
    print("Looking for homologues")
    results = blastq.qblast(
        "blastp", "swissprot", myfasta.seq, alignments=500, word_size=6
    )
    blastRecord = blastparser.read(results)
    #  write hit file:
    with open(f"{gene}_hits.fasta", "w") as f:
        for alignment in blastRecord.alignments:
            f.write(f">{str(myfasta.id)}\n")
            f.write(str(myfasta.seq).strip() + "\n\n")
            for hsp in alignment.hsps:
                f.write(f">{alignment.hit_id}\n")
                f.write(str(hsp.sbjct).replace("-", "") + "\n\n")
    print("Done")


def run_hmm(
    gene_path: Path,
    gene: str,
    res_cutoff: float,
    seqid_cut: float,
    hmmer_home: str,
    cobalt_home: str,
    max_pdbs: int,
):

    with open("master_isoform.txt", "r") as f:
        iso_list = f.read().strip().splitlines()
    isonames = [i.replace(".fasta", "") for i in iso_list]

    with open("usable_isoforms.txt", "w") as f:
        isonames = [
            iso for iso in isonames if Path(iso + ".fasta").stat().st_size != 0
        ]
        for iso in isonames:
            outdir = make_isof_folders(gene_path, iso)
            myfasta = SeqIO.read(
                Path(outdir.parent.absolute(), f"{str(iso)}.fasta"), "fasta"
            )
            start_query(myfasta, gene)
            time.sleep(1)
            check_cobalt = run_cobalt(gene, cobalt_home)
            if not check_cobalt:
                build_hmm(gene, hmmer_home)
                check_output = search_templates(gene, max_pdbs)
                if check_output:
                    print(
                        "We are probably timed out, we'll wait 5 minutes and try again"
                    )
                    time.sleep(300)
                    check_output = search_templates(gene, max_pdbs)
                if not check_output:
                    pdbs, chains = get_templates()
                    toremove = []
                    for pdbs_to_remove_positions_counter, (
                        pdb,
                        chain,
                    ) in enumerate(zip(pdbs, chains)):
                        pdb_exist = download_pdb(pdb)
                        if pdb_exist:
                            check_pdb = manage_pdb(
                                pdb, chain
                            )  # NMR structures (no resolution) will be removed
                        else:
                            check_pdb = False
                        if not check_pdb:
                            toremove.append(pdbs_to_remove_positions_counter)
                    if toremove:
                        for counter, position in enumerate(toremove):
                            pdbs.pop(position - counter)
                            chains.pop(position - counter)
                    for pdb, chain in zip(pdbs, chains):
                        get_fasta(pdb, chain)
                    merge_fastas(pdbs, chains, gene)
                    align_templates_to_msa(gene, cobalt_home)
                    parser_aligned("final_msa.fasta", slash)
                    used_beginnings, used_ends, used_chains, used_pdbs = (
                        run_templates_analysis(res_cutoff, seqid_cut, slash)
                    )
                    if used_beginnings == 42:
                        check_output = True
                os.chdir("..")
                if not check_output:
                    f.write(iso)
                    if iso != iso[-1]:
                        f.write("\n")
            else:
                os.chdir("..")
        else:
            check_output = True
    return check_output
