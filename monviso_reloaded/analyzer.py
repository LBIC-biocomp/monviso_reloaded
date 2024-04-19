import os
import sys
import numpy as np
import torch as pt 
import pandas as pd
from tqdm import tqdm
from glob import glob
from pathlib import Path
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB import PDBParser
from .file_handler import FileHandler


class Analyzer:
    # pylint: disable=import-error
    def __init__(self,pesto_path,output_path,gene_list):
        self.output_path=output_path
        self.gene_list=gene_list
        model_path=os.path.join(pesto_path,"model/save/i_v4_1_2021-09-07_11-21/")
        self.model_filepath=os.path.join(model_path,"model_ckpt.pt")
        if pesto_path not in sys.path:
            sys.path.insert(0, pesto_path)
        if model_path not in sys.path:
            sys.path.insert(0, model_path)
            
        from config import config_model # type: ignore
        from model import Model # type: ignore
        
        self.device=pt.device("cpu")
        self.model = Model(config_model)
        self.model.load_state_dict(pt.load(self.model_filepath, map_location=pt.device("cpu")))
        self.model = self.model.eval().to(self.device)
        
        self._findPDBfiles()
        self.runPestoAnalysis()
        self.runSasaAnalysis()
    

    def _findPDBfiles(self):
        
        from src.dataset import StructuresDataset # type: ignore
        self.pdb_filepaths=[]
        for gene in self.gene_list:
            self.pdb_filepaths+=glob(str(Path(self.output_path,gene,"*","*model","*pdb")))
        
        self.pdb_filepaths=[f for f in self.pdb_filepaths if "pesto" not in f]
        self.dataset=StructuresDataset(self.pdb_filepaths, with_preprocessing=True) # type: ignore
    
    def runPestoAnalysis(self):
        from src.structure import concatenate_chains # type: ignore
        from src.data_encoding import encode_structure, encode_features, extract_topology # type: ignore
        from src.dataset import collate_batch_features # type: ignore
        from src.structure import encode_bfactor,split_by_chain # type: ignore
        from src.structure_io import save_pdb # type: ignore
        
        #module 'numpy' has no attribute 'object'. The aliases was originally deprecated in NumPy 1.20
        #This fixes compatibility with the Numpy versio PeSto was programmed for.
        setattr(np,"object",object) 
        
        results=["Protein", "DNA-RNA", "Lipid", "Ligand", "Ion"]
        with pt.no_grad():
            for subunits, filepath in tqdm(self.dataset):
                # concatenate all chains together
                structure = concatenate_chains(subunits)

                # encode structure and features
                X, M = encode_structure(structure) 
                #q = pt.cat(encode_features(structure), dim=1)
                q = encode_features(structure)[0] 

                # extract topology
                ids_topk, _, _, _, _ = extract_topology(X, 64)

                # pack data and setup sink (IMPORTANT)
                X, ids_topk, q, M = collate_batch_features([[X, ids_topk, q, M]]) 

                # run model
                z = self.model(X.to(self.device), ids_topk.to(self.device), q.to(self.device), M.float().to(self.device))

                # for all predictions
                for i in range(z.shape[1]):
                    # prediction
                    p = pt.sigmoid(z[:,i])

                    # encode result
                    structure = encode_bfactor(structure, p.cpu().numpy())

                    # save results
                    output_filepath = filepath[:-4]+'_pesto_{}.pdb'.format(results[i])
                    save_pdb(split_by_chain(structure), output_filepath)
    
    def runSasaAnalysis(self) -> None:
        p = PDBParser(QUIET=1)
        with FileHandler() as fh:
            for pdb in self.pdb_filepaths:
                
                struct_sasa_path=Path(pdb.replace(".pdb",".struct.sasa"))
                residue_sasa_path=Path(pdb.replace(".pdb",".residue.sasa.csv"))

                struct = p.get_structure("Protein", pdb)
                sr = ShrakeRupley()
                sr.compute(struct, level="S")
                sr.compute(struct, level="R")
                sr.compute(struct, level="A")

                struct_sasa_string=f"Strcuture SASA: {struct.sasa}"
                fh.write_file(struct_sasa_path,struct_sasa_string)
                
                residue_sasa=[]
                backbone_sasa=[]
                sidechain_sasa=[]
                
                for res in struct.get_residues():
                    residue_sasa.append(res.sasa)
                    bb_atoms=[a for a in res.get_atoms() if a.name in ['C','N','CA','O']]
                    sc_atoms=[a for a in res.get_atoms() if a.name not in ['C','N','CA','O']]
                    backbone_sasa.append(sum([a.sasa for a in bb_atoms]))
                    sidechain_sasa.append(sum([a.sasa for a in sc_atoms]))
                    
                residue_sasa_string="Residue sasa, Backbone, Sidechain\n"
                
                for i,r in enumerate(residue_sasa):
                    residue_sasa_string+=f"{r},{backbone_sasa[i]},{sidechain_sasa[i]}\n"
                
                fh.write_file(residue_sasa_path,residue_sasa_string)