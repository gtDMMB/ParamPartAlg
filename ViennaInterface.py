from scorer import RNAStructure
import os
import subprocess
import RNA

class ViennaInterface:
    def __init__(self, filePath:str, transform: bool = False, dangles = 2, noLP = True): 
        self.filePath = filePath
        self.name = os.path.splitext(os.path.basename(filePath))[0]
        self.transform = transform
        self.mfeCalls = 0
        self.suboptCalls = 0
        self.truePmfeCalls = 0 #Calls without lookup
        self.seq = self.sequence_from_fasta()
        self.dangles = dangles
        self.noLP = noLP

        if dangles != 2:
            self.mfe_func = lambda a, b, c, d: self.get_mfe_d1(a, b, c, d)
            try:
                os.mkdir("temp_params")
            except FileExistsError:
                pass
        else:
            self.mfe_func = lambda a, b, c, d: self.get_mfe(a, b, c, d)
            #Setup Vienna
            RNA.params_load_RNA_Turner2004()
            self.md = RNA.md()
            self.md.uniq_ML = 1
            self.md.compute_bpp = 0
            self.md.noLP = noLP
            self.md.dangles = self.dangles

            self.fc = RNA.fold_compound(self.seq, self.md)
            self.mfe_func
    
        #Save computed points
        self.computed = {} #param:sig

    #Calls pmfe with params a, b, c, d
    def vertex_oracle(self, a: int, c: int, b=0, d=1):
        self.mfeCalls += 1
        try:
            return self.computed[(a, b, c, d)][0]
        except KeyError:
            self.truePmfeCalls += 1
            if self.transform:
                a1 = a-(c*3)
                mfe = self.mfe_func(a1, b, c, d)
                sig, struct = mfe
                transformed_sig = self.transform_z(sig)
                self.computed[(a, b, c, d)] = transformed_sig, struct
                return transformed_sig #self.transform_z

            mfe = self.mfe_func(a,b,c,d)
            sig, struct = mfe
            self.computed[(a, b, c, d)] = sig, struct
            return sig

    #Calls subopt with params a, b, c, d
    def subopt_oracle(self, a: int, b: int, c: int, d: int, eng=15):
        self.suboptCalls += 1
        if self.transform:
            a1 = a-(c*3)
            subopt= self.get_subopt(a1,b,c,d)
            transformed = []
            for s in subopt:
                transformed.append(self.transform_z(s))

            return transformed
        return self.get_subopt(a,b,c,d,eng=eng)

#----------------------------------------------------------------------------------------------------
    # Internal Methods
        
    def compute_w(self, sig, params, eng):
        #times 100 to convert from kkals to dckals
        return (round(eng*100) - sum((s*p) for s,p in zip(sig, params[:3])))/params[3]

    def get_mfe_d1(self, a, b, c, d=1):
        full_param_string = lambda a, b, c: f"""## RNAfold parameter file v2.0\n# ML_params\n{b} 0 {a} 0 {c} 0"""
        f_name = os.path.join(f"temp_params", f"{self.name}_tMLparams_LP{"n" if self.noLP else "y"}_d{self.dangles}.par")
        with open(f_name, "w") as f:
            f.write(full_param_string(a, b, c))
        
        command_str = f"RNAfold --infile {self.filePath} --paramFile {f_name} -d {self.dangles}"
        if self.noLP:
            command_str = command_str + " --noLP"
        
        command = command_str.split()

        run_command = subprocess.run(command, stdout=subprocess.PIPE, encoding='UTF-8')
        seq, structure, en = run_command.stdout.replace(" ( ", " (").split()

        en = eval(en)

        sig_no_w = tuple(RNAStructure(structure).score_structure())
        sig = sig_no_w + (self.compute_w(sig_no_w, (a,b,c,1), en),)
        return sig, structure

    def change_params(self, a, b, c):
        #String b 0 a 0 c 0 found from /ViennaRNA-2.5.1/misc/rna_turner1999.par (search for "ML_params")
        string = f"## RNAfold parameter file v2.0\n# ML_params\n{b} 0 {a} 0 {c} 0 \n"

        RNA.params_load_from_string(string)
        p = RNA.param(self.md)
        exp_p = RNA.exp_param(self.md)

        # substitute energy parameters within the fold_compound
        # with the 'new' global settings
        self.fc.params_subst(p)
        self.fc.exp_params_subst(exp_p)

    #NEED TO TEST THIS METHOD RIGOROUSLY
    def get_mfe(self, a, b, c, d=1):
        self.change_params(a, b, c)
        # compute MFE
        (structure, en) = self.fc.mfe()
        sig_no_w = tuple(RNAStructure(structure).score_structure())

        sig = sig_no_w + (self.compute_w(sig_no_w, (a,b,c,d), en),)
        return sig, structure
    
    def sequence_from_fasta(self):
        with open(self.filePath) as f:
            seq = f.readline()
            while seq[0] == ">":
                seq = f.readline()
            return seq.strip("\n")

    def get_subopt(self, a, b, c, d=1, eng=5):
        a = int(round(a, 0))
        b = int(round(b, 0))
        c = int(round(c, 0))

        self.change_params(a, b, c)

        # compute MFE
        sigs_no_w = ((tuple(RNAStructure(s.structure).score_structure()), s.energy) for s in self.fc.subopt(eng))
        subopt = [s[:3] + (self.compute_w(s, (a,b,c,d),en),) for s, en in sigs_no_w]

        return subopt
    
    def transform_z(self, point):
        x,y,z,w = point
        return (x, y, z - (3*x), w)