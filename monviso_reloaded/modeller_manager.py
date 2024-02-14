from pathlib import Path
from .file_handler import FileHandler

class Modeller_manager:
    def __init__(self,isoform,mutation):
        self.isoform=isoform
        self.mutation=mutation
        self.sequence_to_model=self.isoform.aligned_sequence[:]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass
    
    def run(self):
        
        print(f"Modelling {self.isoform.gene_name} {self.isoform.isoform_name} {self.mutation}")
        if self.mutation!="WT":
            success_mutation=self._mutate_reside(self.mutation)
            if not success_mutation:
                print(self.sequence_to_model)
                raise(RuntimeError(f"Could not apply mutation!"))
        self.write_alignment()
        self.write_script()
            


    def write_script(self):
        alignment_name="modeller_input_"+"".join(self.mutation)+".dat"
        output_name="modeller_output_"+"".join(self.mutation)+".dat"
        script_name="run_modeller_"+"".join(self.mutation)+".py"
        template_names=[t.pdb_name+"_"+t.pdb_chain+"_clean" for t in self.isoform.templates]
        content="""from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.parallel import *
import os

log.verbose()
j=job()
cores_number = 8
thread = 0
while thread < cores_number:
    j.append(local_slave())
    thread += 1
env = environ()
env.io.atom_files_directory = """+ str(Path(self.isoform.out_path,"templates"))+"""
env.io.hetatm = True
a = automodel(env,
    alnfile =\""""+alignment_name+"""\",
    knowns = ("""+ str(template_names)+"""),
    sequence = \""""+self.isoform.gene_name+"""\",
    assess_methods=(assess.DOPE, assess.GA341))
a.starting_model= 1
a.ending_model  = 1
a.use_parallel_job(j)
a.make()
ok_models = filter(lambda x: x['failure'] is None, a.outputs)
toscore = 'DOPE score'
ok_models = sorted(ok_models, key=lambda k: k[toscore])
models = [m for m in ok_models[0:10]]
myout = open(\""""+output_name+"""\", "w")
for m in models:
        myout.write(str(m['name']) + " (DOPE SCORE: %.3f)" % (m[toscore]))
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')
mdl = complete_pdb(env, m['name'])
s = selection(mdl)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=\""""+self.isoform.gene_name+"""\", normalize_profile=True, smoothing_window=15)""" 
       
        with FileHandler() as fh:
            fh.write_file(Path(self.isoform.out_path,script_name),content)
            
    def _mutate_reside(self,mutation) -> bool:
        """Take a mutation in the format [1 letter amino acid, residue number, 1 lett. amino acid]
        and apply it to the aligned sequence to model, taking into account all the "-"'s.
        
        Returns True, if succesful, else False
        """
        i=int(mutation[1])-1
        while i<len(self.sequence_to_model):
            actual_residue_index=i-self.sequence_to_model[:i].count("-")
            if  (actual_residue_index+1==int(mutation[1])) and self.sequence_to_model[i]==mutation[0]:
                self.sequence_to_model=self.sequence_to_model[:i]+mutation[2]+self.sequence_to_model[i+1:]
                return True
            i+=1
        
        print(f"Could not apply mutation {mutation}!")
        return False
    
    def write_alignment(self):
        
        alignment_name="modeller_input_"+"".join(self.mutation)+".dat"
        output_path=Path(self.isoform.out_path,alignment_name)
        # Create an object storing all names and aligned sequences
        sequences=[[self.isoform.gene_name,self.sequence_to_model]]
        for template in self.isoform.templates:
            sequences.append([template.pdb_name+"_"+template.pdb_chain,template.aligned_sequence])
            
        ## TO DO here: add chain breaks in place of long seqs with no coverage
        # <--------   
        #
            
        # Start writing the content string to be printed in the file
        content=""
        content+=">P1;"+sequences[0][0]+"\n"
        content+="sequence:"+sequences[0][0]+"\n"
        content+=sequences[0][1]+"*\n"
        
        # Add templates
        for template_sequence in sequences[1:]:
            content+=">P1;"+template_sequence[0]+"_clean"+"\n"
            content+="structureX:"+template_sequence[0]+"_clean"+"\n"
            content+=template_sequence[1]+"*\n"
            
        with FileHandler() as fh:
            fh.write_file(output_path,content)