from pathlib import Path
from .file_handler import FileHandler
import subprocess
import pandas as pd

class DockableStructure:
    def __init__(self, gene, path):
        self.gene=gene
        self.path=path
        self.name=path.name
        self._load_analysis_files()
    
    def _load_analysis_files(self):
        self.pestoprotein=Path(str(self.path).replace(".pdb","_pesto_Protein.pdb"))
        self.residuesasa=Path(str(self.path).replace(".pdb",".residue.sasa.csv"))
        self.residuedepth=Path(str(self.path).replace(".pdb",".residue.depth.csv"))
        with FileHandler() as fh:
            check1=fh.check_existence(self.pestoprotein)
            check2=fh.check_existence(self.residuesasa)
            check3=fh.check_existence(self.residuedepth)
        if sum([check1,check2,check3])==3:
            pass
        else:
            raise FileNotFoundError(f"Pesto and residues analysis was not completed for gene {self.gene}.")
        
        self.residuesasa_db=pd.read_csv(self.residuesasa)
        
    def change_path(self,path_to: Path):
        with FileHandler() as fh:
            fh.copy_file(self.path,Path(path_to,self.name))
            self.path=Path(path_to,self.name)
    
    def write_sasa_residues(self, output_path):
        pass
            
class HaddockManager:
    def __init__(self, protein1: DockableStructure,protein2: DockableStructure,output:Path):
        self.protein1=protein1
        self.protein2=protein2
        self.output=output
        self.default_config="""# directory name of the run
run_dir = "$runname$"

# compute mode
mode = "local"


# Self contained rundir (to avoid problems with long filename paths)
self_contained = true

# molecules to be docked
molecules =  [ $moleculesname$ ]

[topoaa]

[rigidbody]
# CDR to surface ambig restraints
ambig_fname = "$ambigfilename$"
# Restraints to keep the antibody chains together
unambig_fname = "$unambigfilename$"
# Number of models to generate
sampling = 100

[seletopclusts]
## select the best 10 models of each cluster
top_models = 10
"""
    def copy_structures(self):
        """Copy the pdb structure of the two proteins in the docking folder.
        Change the chain of the second protein from "A" to "B".
        """
        self.protein1.change_path(self.output)
        self.protein2.change_path(self.output)
        
    def find_most_exposed_residues(self):
        pass

class DockingManager:
    def __init__(self,output_path,gene_list,haddock_home,hdocklite_home,megadock_home):
        self.output_path=Path(output_path)  
        self.gene_list=gene_list
        self.haddock_home=haddock_home
        self.hdocklite_home=hdocklite_home
        self.megadock_home=megadock_home
        
    def run(self):
        self.load_structures()
        self.make_folders()
        self.run_megadock()
        self.run_hdocklite()
        
    def load_structures(self):
        """For each couple of genes in self.gene_list,
        search the corresponding folder for modelled structures.
        """
        self.coupled_structure_lists=[]
        
        for gene_couple in self.gene_list:
            structs=[[],[]]
            pattern = f"{gene_couple[0]}/isoform*/{gene_couple[0]}*model/*.pdb*"
            structures=list(self.output_path.glob(pattern))
            structs[0] = [DockableStructure(gene_couple[0],s) for s in structures if "pesto" not in str(s)]
            pattern = f"{gene_couple[1]}/isoform*/{gene_couple[1]}*model/*.pdb*"
            structures=list(self.output_path.glob(pattern))
            structs[1] = [DockableStructure(gene_couple[1],s) for s in structures if "pesto" not in str(s)]
            
            self.coupled_structure_lists.append(structs)

    def make_folders(self):
        with FileHandler() as fh:
            fh.create_directory(Path(self.output_path,"Docked"))
            for gene_couple in self.gene_list:
                fh.create_directory(Path(self.output_path,"Docked","-".join(gene_couple)))
                fh.create_directory(Path(self.output_path,"Docked","-".join(gene_couple),"MEGADOCK"))
                fh.create_directory(Path(self.output_path,"Docked","-".join(gene_couple),"HDOCKLITE"))
                fh.create_directory(Path(self.output_path,"Docked","-".join(gene_couple),"HADDOCK"))
    
    
    def run_megadock(self,n_exported_structs=100):
        for couple in self.coupled_structure_lists:
            for p1 in couple[0]:
                for p2 in couple[1]:
                    directory_name=p1.gene+"-"+p2.gene
                    file_name=p1.name+"-"+p2.name
                    file_name=file_name.replace(".pdb","")+".out"

                    output=Path(self.output_path,"Docked",directory_name,"MEGADOCK",file_name)
                    with FileHandler() as fh:
                        if not fh.check_existence(output):
                            command = f"{str(Path(self.megadock_home,'megadock'))} -R {str(p1.path)} -L {p2.path} -o {str(output)}"
                            subprocess.run(
                                command, shell=True, universal_newlines=True, check=True
                            )
                    
                        for i in range(n_exported_structs):
                            exported_pdb=file_name.replace(".out",f".{i}.pdb")
                            exported_pdb_path=Path(self.output_path,"Docked",directory_name,"MEGADOCK",exported_pdb)

                            if not fh.check_existence(exported_pdb_path):
                                command = f"{str(Path(self.megadock_home,'decoygen'))} {str(exported_pdb_path)} {p2.path} {str(output)} {i+1}"
                                subprocess.run(
                                command, shell=True, universal_newlines=True, check=True
                                )
    
    def run_hdocklite(self,n_exported_structs=100):
        for couple in self.coupled_structure_lists:
            for p1 in couple[0]:
                for p2 in couple[1]:
                    directory_name=p1.gene+"-"+p2.gene
                    file_name=p1.name+"-"+p2.name
                    file_name=file_name.replace(".pdb","")+".out"

                    output=Path(self.output_path,"Docked",directory_name,"HDOCKLITE",file_name)
                    with FileHandler() as fh:
                        if not fh.check_existence(output):
                            command = f"{str(Path(self.hdocklite_home,'hdock'))} {str(p1.path)} {p2.path} -out {str(output)}"
                            subprocess.run(
                                command, shell=True, universal_newlines=True, check=True
                            )
                    
                        
                        exported_pdb=file_name.replace(".out",f".top{n_exported_structs}.pdb")
                        exported_pdb_path=Path(self.output_path,"Docked",directory_name,"HDOCKLITE",exported_pdb)

                        if not fh.check_existence(exported_pdb_path):
                            command = f"{str(Path(self.hdocklite_home,'creapl'))} {str(output)} {str(exported_pdb_path)} -nmax {n_exported_structs} -complex -models"
                            subprocess.run(
                            command, shell=True, universal_newlines=True, check=True
                            )
    
    def run_haddock(self):
        hm= HaddockManager()
        for couple in self.coupled_structure_lists:
            for p1 in couple[0]:
                for p2 in couple[1]:
                    directory_name=p1.gene+"-"+p2.gene
                    file_name=p1.name+"-"+p2.name
                    file_name=file_name.replace(".pdb","")+"_run"

                    output=Path(self.output_path,"Docked",directory_name,"HADDOCK",file_name)
                    with FileHandler() as fh:
                        if not fh.check_existence(output):
                            configfile=hm.createConfig()
                            hm.run(configfile)
                    
