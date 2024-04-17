import os
import sys
import h5py
import json
import numpy as np
import torch as pt 
import pandas as pd
from tqdm import tqdm
from glob import glob
from pathlib import Path

class PestoAnalyzer:
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
        self.run()
    

    def _findPDBfiles(self):
        
        from src.dataset import StructuresDataset # type: ignore
        pdb_filepaths=[]
        for gene in self.gene_list:
            pdb_filepaths+=glob(str(Path(self.output_path,gene,"*","*model","*pdb")))
        self.dataset=StructuresDataset(pdb_filepaths, with_preprocessing=True) # type: ignore
    
    def run(self):
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