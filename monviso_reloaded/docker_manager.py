import subprocess
import pandas as pd
import os
import numpy as np
from pathlib import Path
from scipy.spatial.transform import Rotation

from .file_handler import FileHandler
from .PDB_manager import PDB_manager

class DockableStructure:
    def __init__(self, gene, path):
        self.gene=gene
        self.path=path
        self.name=path.name
        self.cutoff=90
        self._load_analysis_files()
    
    def set_cutoff(self,c):
        self.cutoff=float(c)
    
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
                
    def change_path(self,path_to: Path,new_chain_name="A"):
        with PDB_manager() as pm:
            pm.change_chain_and_save(str(self.path),new_chain_name,str(path_to))
            self.path=Path(path_to)
    
    def write_sasa_residues(self, output_path):
        "Writes residues in the top10% percentile for absolute SASA."
        self.residuesasa_db=pd.read_csv(self.residuesasa)
        cutoff=(100-self.cutoff)/100
        top10threshold=self.residuesasa_db[' Residue sasa'].quantile(cutoff)
        topresidues=list(self.residuesasa_db[self.residuesasa_db[' Residue sasa']>top10threshold]["Residue number"])
        output_string=",".join([str(x) for x in topresidues])+"\n"
        with FileHandler() as fh:
            fh.write_file(output_path,output_string)
    
    def write_pesto_residues(self, output_path):
        "Writes residues in the top10% percentile for PeSTo predicted score.."
        with FileHandler() as fh:
            pestofile=fh.read_file(self.pestoprotein).split("\n")
            CAs=[line for line in pestofile if "CA" in line]
            resnums=[line[22:26] for line in CAs]
            scores=[line[62:66] for line in CAs]
            df=pd.DataFrame(zip(resnums, scores),columns=["Residue number","PeSTo score"])
            threshold=df["PeSTo score"].map(float).quantile(0.9)
            selection=df[df["PeSTo score"].map(float)>threshold]
            output_string=",".join(list(selection["Residue number"]))
            fh.write_file(output_path,output_string)

class ClusteringStructure:
    def __init__(self,path, software):
        self.path=path
        self.software=software
        self.ca_positions=[]
        self.chain_a_count=0
        self.chain_b_count=0
        self._load_CA_atoms()

    
    def _load_CA_atoms(self):
        ca_positions = []

        with FileHandler() as fh:
            lines = fh.read_file(self.path).split('\n')

        for line in lines:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                ca_positions.append([x, y, z])

                chain_id = line[21].strip()
                if chain_id == 'A':
                    self.chain_a_count += 1
                elif chain_id == 'B':
                    self.chain_b_count += 1

        self.ca_positions = np.array(ca_positions)
    
    def fit_on_reference(self,reference):
        if self.chain_a_count>self.chain_b_count:
            ref_ca_positions=reference.ca_positions[:self.chain_a_count]
            fit_ca_positions=self.ca_positions[:self.chain_a_count]
        else:
            ref_ca_positions=reference.ca_positions[self.chain_b_count:]
            fit_ca_positions=self.ca_positions[self.chain_b_count:]
            
        fit_ca_positions-ref_ca_positions.mean(axis=0)
        self.ca_positions-ref_ca_positions.mean(axis=0)

        rotation_fit=Rotation.align_vectors(ref_ca_positions,fit_ca_positions)
        self.ca_positions=rotation_fit[0].apply(self.ca_positions)
    
    def calculate_smaller_chain_RMSD(self,reference):
        if self.chain_a_count<self.chain_b_count:
            ref_ca_positions=reference.ca_positions[:self.chain_a_count]
            fit_ca_positions=self.ca_positions[:self.chain_a_count]
        else:
            ref_ca_positions=reference.ca_positions[self.chain_b_count:]
            fit_ca_positions=self.ca_positions[self.chain_b_count:]
            
        distances = np.linalg.norm(ref_ca_positions - fit_ca_positions, axis=1)
        distances_sqrd= distances**2
        rmsd= np.sqrt( distances_sqrd.mean())
        return rmsd
        
        
class HaddockManager:
    def __init__(self, protein1: DockableStructure,protein2: DockableStructure,output:Path, haddock_selection: str, haddock_cutoff: float):
        self.protein1=protein1
        self.protein1.set_cutoff(haddock_cutoff)
        self.protein2=protein2
        self.protein2.set_cutoff(haddock_cutoff)
        self.output=output
        self.haddock_selection=haddock_selection
        self.haddock_cutoff=haddock_cutoff
        self.default_config="""
# directory name of the run
run_dir = "run"

# compute mode
mode = "local"
#  5 nodes x 50 tasks = 250
ncores = 250

# Self contained rundir
self_contained = true

# Post-processing to generate statistics and plots
postprocess = true

# Cleaninog
clean = true

# molecules to be docked
molecules =  [
	"$molecule1",
    "$molecule2"
    ]

# ====================================================================
# Parameters for each stage are defined below, prefer full paths
# ====================================================================
[topoaa]

[rigidbody]
# paratope to surface ambig restraints
ambig_fname = "$ambig_restraints"
# Turn off ramdom removal of restraints
randremoval = false
# Number of models to generate
sampling = 1000

"""

    def copy_structures(self):
        """Copy the pdb structure of the two proteins in the docking folder.
        Change the chain of the second protein from "A" to "B".
        """
        self.protein1.change_path(self.output)
        self.protein2.change_path(self.output)
        
    
    def generate_tbl(self, residue_path1, residue_path2, output_path, chain1="A",chain2="B"):
        with FileHandler() as fh:
            res1=fh.read_file(residue_path1).split(",")
            res2=fh.read_file(residue_path2).split(",")
            
            output_string=""
            
            r1selections=[f"(resid {r1} and segid {chain1})" for r1 in res1]
            r2selections=[f"(resid {r2} and segid {chain2})" for r2 in res2]
            
            for i,r1 in enumerate(r1selections):
                output_string+="\nassign "+r1+"\n(\n"
                output_string+="\nor\n".join(r2selections)
                output_string+=")  2.0 2.0 0.0"
            
            fh.write_file(output_path,output_string)
            
    def generate_config(self, protein1: DockableStructure, protein2: DockableStructure, tbl_path: Path, output_path: Path):
        output_string=self.default_config
        output_string=output_string.replace("$molecule1",str(protein1.path.absolute()))
        output_string=output_string.replace("$molecule2",str(protein2.path.absolute()))
        output_string=output_string.replace("$ambig_restraints",tbl_path.name)

        with FileHandler() as fh:
            fh.write_file(output_path,output_string)
            
    def run(self, config_path: Path):
        cwd=os.getcwd()
        os.chdir(config_path.parent)
        
        command = f"haddock3 {config_path.name}"
        subprocess.run(command, shell=True, universal_newlines=True, check=True)
        os.chdir(cwd)
        

class DockingManager:
    def __init__(self,output_path,gene_list,haddock_home,haddock_selection,haddock_cutoff,hdocklite_home,megadock_home):
        self.output_path=Path(output_path)  
        self.gene_list=gene_list
        self.haddock_home=haddock_home
        self.haddock_selection=haddock_selection
        self.haddock_cutoff=haddock_cutoff
        self.hdocklite_home=hdocklite_home
        self.megadock_home=megadock_home
        
    def run(self):
        self.load_structures()
        self.make_folders()
        self.run_megadock()
        self.run_hdocklite()
        self.run_haddock()
        self.cluster_structures()
        #self.run_final()
        
    def load_structures(self):
        """For each couple of genes in self.gene_list,
        search the corresponding folder for modelled structures.
        """
        self.coupled_structure_lists=[]
        
        for gene_couple in self.gene_list:
            structs=[[],[]]
            pattern = f"{gene_couple[0]}/*/{gene_couple[0]}*model/*.pdb*"
            structures=list(self.output_path.glob(pattern))
            structs[0] = [DockableStructure(gene_couple[0],s) for s in structures if "pesto" not in str(s)]
            pattern = f"{gene_couple[1]}/*/{gene_couple[1]}*model/*.pdb*"
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
                fh.create_directory(Path(self.output_path,"Docked","-".join(gene_couple),"FINAL"))
    
    
    def run_megadock(self,n_exported_structs=100):
        for couple in self.coupled_structure_lists:
            for p1 in couple[0]:
                for p2 in couple[1]:
                    directory_name=p1.gene+"-"+p2.gene
                    file_name=p1.name+"-"+p2.name
                    file_name=file_name.replace(".pdb","")

                    output=Path(self.output_path,"Docked",directory_name,"MEGADOCK",file_name,file_name+".out")
                    with FileHandler() as fh:
                        fh.create_directory(Path(self.output_path,"Docked",directory_name,"MEGADOCK",file_name))
                        if not fh.check_existence(output):
                            
                            command = f"{str(Path(self.megadock_home,'megadock'))} -R {str(p1.path)} -L {p2.path} -o {str(output)}"
                            subprocess.run(
                                command, shell=True, universal_newlines=True, check=True
                            )

                            for i in range(n_exported_structs):
                                p2.change_path(p2.path,"B")
                                exported_pdb=file_name+f".{i}.pdb"
                                exported_pdb_path=Path(self.output_path,"Docked",directory_name,"MEGADOCK",file_name,exported_pdb)

                                if not fh.check_existence(exported_pdb_path):
                                    command = f"{str(Path(self.megadock_home,'decoygen'))} {str(exported_pdb_path)} {p2.path} {str(output)} {i+1}"
                                    subprocess.run(
                                    command, shell=True, universal_newlines=True, check=True
                                    )
                                p2.change_path(p2.path,"A")

                                megadock_generated_pdb=fh.read_file(exported_pdb_path)
                                protein_partner= fh.read_file(p1.path)
                                content=protein_partner+megadock_generated_pdb

                                fh.write_file(exported_pdb_path, content)
                        else:
                            print(f"Skipping MEGADOCK job. Output {str(output)} already exists.")
    
    def run_hdocklite(self,n_exported_structs=100):
        for couple in self.coupled_structure_lists:
            for p1 in couple[0]:
                for p2 in couple[1]:
                    directory_name=p1.gene+"-"+p2.gene
                    file_name=p1.name+"-"+p2.name
                    file_name=file_name.replace(".pdb","")

                    tmp_output=Path("./","Hdock.out")
                    output=Path(self.output_path,"Docked",directory_name,"HDOCKLITE",file_name,file_name+".out")
                    with FileHandler() as fh:
                        fh.create_directory(Path(self.output_path,"Docked",directory_name,"HDOCKLITE",file_name))
                        if not fh.check_existence(output):
                            p2.change_path(p2.path,"B")
                            command = f"{str(Path(self.hdocklite_home,'hdock'))} {str(p1.path)} {p2.path} -out {str(tmp_output)}"
                            subprocess.run(
                                command, shell=True, universal_newlines=True, check=True
                            )
                            p2.change_path(p2.path,"A") #Revert change
                            fh.move_file(tmp_output,output)
                        else:
                            print(f"Skipping HDOCKlite job. Output {str(output)} already exists.")
                        
                        exported_pdb=[file_name+f"_{ID+1}.pdb" for ID in range(n_exported_structs)]
                        exported_pdb_path=[Path(self.output_path,"Docked",directory_name,"HDOCKLITE",file_name,file) for file in exported_pdb]
                        
                        command = f"{str(Path(self.hdocklite_home,'createpl'))} {str(output)} model.pdb -nmax {n_exported_structs} -complex -models"
                        subprocess.run(
                            command, shell=True, universal_newlines=True, check=True
                            )

                        for modelID in range(n_exported_structs):
                                fh.move_file(Path(Path.cwd(),f"model_{modelID+1}.pdb"),exported_pdb_path[modelID])
                        

    def run_haddock(self):
        for couple in self.coupled_structure_lists:
            for p1 in couple[0]:
                for p2 in couple[1]:
                    directory_name=p1.gene+"-"+p2.gene
                    file_name=p1.name+"-"+p2.name
                    file_name=file_name.replace(".pdb","")+"_run"
                    
                    output=Path(self.output_path,"Docked",directory_name,"HADDOCK",file_name)

                    
                    sasa_res1=Path(output,p1.name+"_sasa_residues.txt")
                    pesto_res1=Path(output,p1.name+"_pesto_residues.txt")
                    sasa_res2=Path(output,p2.name+"_sasa_residues.txt")
                    pesto_res2=Path(output,p2.name+"_pesto_residues.txt")
                    
                    sasa_ambig=Path(output,"sasa_ambig.tbl")
                    pesto_ambig=Path(output,"pesto_ambig.tbl")
                    config_path=Path(output,"docking.cfg")

                    with FileHandler() as fh:
                        if not fh.check_existence(output):
                            fh.create_directory(output)
                            hm= HaddockManager(p1,p2,output,self.haddock_selection, self.haddock_cutoff)
                            #configfile=hm.createConfig()
                            p1.write_sasa_residues(sasa_res1)
                            p1.write_pesto_residues(pesto_res1)
                            p2.write_sasa_residues(sasa_res2)
                            p2.write_pesto_residues(pesto_res2)
                            #p1.change_path(Path(output,p1.name),"A")
                            p2.change_path(p2.path,"B")
                            
                            hm.generate_tbl(sasa_res1,sasa_res2,sasa_ambig)
                            hm.generate_tbl(pesto_res1,pesto_res2,pesto_ambig)
                            if self.haddock_selection=="pesto":
                                tbl_path=pesto_ambig
                            elif self.haddock_selection=="sasa":
                                tbl_path=sasa_ambig
                            else:
                                p2.change_path(p2.path,"A") #Revert chain name change
                                raise ValueError("The haddock selection type was not defined correctly. Use \"pesto\" or \"sasa\", without quotation marks.")
                            hm.generate_config(p1,p2,tbl_path,config_path)
                            hm.run(config_path)
                        
                        else:
                            print(f"Skipping Haddock job. Folder {str(output)} already exists.")
                            
                        p2.change_path(p2.path,"A") #Revert chain name change
                        
    
    def cluster_distance_matrix(self,distance_matrix):
        """
        Clustering the structures based on a distance matrix.

        This method applies a clustering algorithm similar to GROMOS for 
        Molecular Dynamics frames, or Buting clustering for molecules,
        to the given distance matrix to group structures into clusters based
        on a cutoff distance. The cutoff distance is optimized automatically
        to have between 70% and 85% of the structures in the first three clusters.

        Args:
        distance_matrix (np.ndarray): Distance matrix between structures.

        Returns:
        list: List of clusters, where each element is a tuple containing the
              centroid index and the list of indices of the structures in the cluster.
        """
        cutoff_max=distance_matrix.max()
        cutoff_min=distance_matrix.min()
        cutoff=(cutoff_max+cutoff_min)/2
        top3cluster= 0.0 #percent of structures
        while top3cluster<0.7 or top3cluster>0.85:
            print(f"Clustering with cutoff {round(cutoff,3)} Å...")
            clusters=[]
            processing_distance_matrix=distance_matrix.copy()
            entriesID=np.array([n for n in range(distance_matrix.shape[0])])
            print("Remaining structures...",end="")
            while processing_distance_matrix.shape[0]>0:
                    print(f"{processing_distance_matrix.shape[0]}",end="...")
                    distances_within_cutoff = processing_distance_matrix <= cutoff
                    new_centroid = np.argmax(distances_within_cutoff.sum(axis=1))
                    centroid_original_id = entriesID[new_centroid]
                    selection = distances_within_cutoff[new_centroid].astype(bool)
                    new_cluster = entriesID[selection] 
                    # Remove rows and columns for selected items
                    mask = ~selection
                    entriesID = entriesID[mask]
                    processing_distance_matrix = processing_distance_matrix[mask][:, mask]
                    clusters.append([centroid_original_id, new_cluster])
            print("")
            if len(clusters) < 3:
                cutoff_max = cutoff
            else:
                top3cluster = sum([len(c[1]) for c in clusters[:3]]) / distance_matrix.shape[0]
                if 0.7 <= top3cluster <= 0.85:
                    break
                elif top3cluster > 0.85:
                    cutoff_max = cutoff
                else:
                    cutoff_min = cutoff

            cutoff = (cutoff_min + cutoff_max) / 2
            print("New attempt with different cutoff.")
        
        print(f"Top 10 cluster size: "+" ".join([str(len(c[1])) for c in clusters[:10]]))
        print(f"Top 3 clusters coverage: {round(top3cluster*100,2)}%")
        return clusters
                
    def cluster_structures(self):        
        for couple in self.coupled_structure_lists:
            for p1 in couple[0]:
                for p2 in couple[1]:
                    dirs=[]
                    directory_name=p1.gene+"-"+p2.gene
                    file_name=p1.name+"-"+p2.name
                    file_name=file_name.replace(".pdb","")
                    dirs.append(Path(self.output_path,"Docked",directory_name,"HADDOCK",file_name+"_run","run","1_rigidbody"))
                    dirs.append(Path(self.output_path,"Docked",directory_name,"HDOCKLITE",file_name))
                    dirs.append(Path(self.output_path,"Docked",directory_name,"MEGADOCK",file_name))


                    #Load all .pdb and .pdb.gz files as MDanalysis universes.
                    structures_for_clustering=[] 
                    software_names=["HADDOCK","HDOCKLITE","MEGADOCK"]
                    for d,dir in enumerate(dirs):
                        structures=list(dir.glob("*pdb*"))
                        structures.sort()
                        for s,struct in enumerate(structures):
                            if struct.suffixes == ['.pdb','.gz']:
                                with gzip.open(str(struct), 'rb') as f_in:
                                    with open(str(struct).replace('.pdb.gz','.pdb'), 'wb') as f_out:
                                        shutil.copyfileobj(f_in, f_out)
                                struct = Path(str(struct).replace('.pdb.gz','.pdb'))
                            structures_for_clustering.append( ClusteringStructure(struct,software_names[d]) )

                    
                    #Fix chain A and chain B atom counts for HOCKlite molecules
                    #Copy the number from haddock structures to all other structures
                    for struct in structures_for_clustering:
                        struct.chain_a_count=structures_for_clustering[0].chain_a_count
                        struct.chain_b_count=structures_for_clustering[0].chain_b_count
                        
                    #And then fit all structures on the first structure of Haddock
                    
                    for struct in structures_for_clustering:
                        struct.fit_on_reference(structures_for_clustering[0])

                    #Create distance matrix
                    distance_matrix= np.zeros((len(structures_for_clustering),len(structures_for_clustering)))
                    for s1, struct1 in enumerate(structures_for_clustering):
                        for s2, struct2 in enumerate(structures_for_clustering[s1:]):
                            rmsd=struct1.calculate_smaller_chain_RMSD(struct2)
                            distance_matrix[s1,s1+s2]=rmsd
                            distance_matrix[s1+s2,s1]=rmsd
        
                    self.cluster_distance_matrix(distance_matrix)