from pathlib import Path
from .file_handler import FileHandler
import subprocess

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

class Docker:
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
    
    
    def run_megadock(self):
        for couple in self.coupled_structure_lists:
            for p1 in couple[0]:
                for p2 in couple[1]:
                    directory_name=p1.gene+"-"+p2.gene
                    file_name=p1.name+"-"+p2.name
                    output=Path(self.output_path,"Docked",directory_name,"MEGADOCK",file_name)
                    with FileHandler() as fh:
                        if not fh.check_existence(output):
                            command = f"{str(Path(self.megadock_home,'megadock'))} -R {str(p1.path)} -L {p2.path} -o {str(output)}"
                            subprocess.run(
                                command, shell=True, universal_newlines=True, check=True
                            )