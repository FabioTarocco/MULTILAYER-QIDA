from utils.input import Input
import numpy as np
from qiskit.circuit import QuantumCircuit, Parameter, ParameterVector, Gate
from qiskit.circuit.parametervector import ParameterVectorElement
from qiskit.circuit.library import TwoLocal
from pprint import pprint
import networkx as nx
import warnings


#from qiskit_nature.circuit.library.initial_states import HartreeFock

#_ANSATZ_LIST = ['two_local','CGR','SU_N','SU_N_chem','SSMI','Cphase','Cphase_H','Cphase_fixed_phase','Cphase_fixed_theta','Cphase_Ry','Cphase_Ry_fixed_phase','Cphase_Ry_fixed_theta','toffoli','symetric_map','symetric_group_map','symetric_particle_conserving','particle_conserving_ansatz','symetric_CNOTs','sym_CNOTs_part_conserving','CGR_rz','CGR_rx']

_PARAM_ENT_GATE_LIST = ['crx', 'cry', 'crz', 'cp', 'rzx','su4']
_ENT_GATE_LIST = ['cx', 'cy', 'cz', 'dcx', 'ch', 'swap', 'iswap'] #, 'cs', 'csd', 'csx', 'rxx', 'ryy', 'rzz']
_ANSATZ_ENTANGLEMENT_MAP_LIST = ['full', 'linear', 'reverse_linear', 'pairwise', 'circular', 'sca']

class Ansatz:
    """
    Generate a Ansatz object.

    Attributes:
        input: (Input)  Input object that contains the initial settings.
        n_qubits: (int | None) Size of the Quantum Circuit that have to be created. It can be passed as external input or given in the input.nml file.
        starting_occupations: (int | list) It corresponds to num_orbitals and it used to initialize the initial state.
        reduce_parametes: (bool) If True, the ansatz will only uses one parameter for the qubits that do not appear in the entangler_map.
        initial_state: (str | None) A string that indicates if the initial state is empty or is Hartree-Fock initial state.
        var_form: (qiskit.QuantumCircuit) The circuit created.
        rotation_blocks: (list) List of parametrized rotation that have to be added at each layer of the Parametrized Quantu Circuit.
        type: (str) This is the type of ansatz that we want to create:  heuristic (completely costumable), two_local (qiskit builder).
        entanglement_blocks: (str) Entanglement block that has to be used following the entangler_map.
        entangler_map: (list) List of pairs of qubit that is used to entangle qubits in the circuit. The structure [[1st layer entangle map[pair of qubits]], [2nd layer entangle map[...]], ...]
        additional_entangler_map: (list) Additional entangler map that can be added on each layer.
        parameters: (list) List of parameters that is associated with the parametrized gates used in the circuit.
        depth: (int) Number of layer of the ansatz.
        string: (str) Counter used to correctly number the parameters.
        
    Methods:
        __init__(input: Input, n_qubits: int|None, starting_occupations: int|list)
        build_ansatz():
        two_local():
        heuristic():
    """
    def __init__(self, input: Input, n_qubits:int , bitString:(bool,int)=(False,0), reduce_parameters:bool = False, initial_state = None, ansatz_type = None, entangler_map = None, parameters = None, depth = None, starting_occupations: int|list = 0):
        
        
        if input is None:
            raise Exception("Input is required to construct an Ansatz object!")
        
        self.input = input
        if self.input.ansatz_n_qubit is None and n_qubits is None:
            raise Exception("Number of qubits is required. \n You can insert it in the input file or as external variable.")
        
        if self.input.ansatz_n_qubit is not None:
            self.n_qubits = self.input.ansatz_n_qubit

        if n_qubits is not None:
            self.n_qubits = n_qubits

        if initial_state is None:
            self.initial_state = self.input.ansatz_initial_state
        else:
            self.initial_state = initial_state
        if ansatz_type is None:
            self.type = self.input.ansatz_type
        else:
            self.type = ansatz_type
        if entangler_map is None:
            self.entangler_map = self.input.ansatz_entangler_map
        else:
            self.entangler_map = entangler_map
        self.rotation_blocks = self.input.ansatz_rotation_blocks
        
        self.entanglement_blocks = self.input.ansatz_entanglement_blocks
        self.additional_entangler_map = self.input.ansatz_additional_entangler_map
        if parameters is None:
            self.parameters = self.input.ansatz_parameters
        else:
            self.parameters = parameters
        if depth is None:
            self.depth = self.input.ansatz_depth
        else:
            self.depth = depth
        self.starting_occupations = starting_occupations
        self.reduce_parameters = reduce_parameters
        
        if self.reduce_parameters:
            self._entangler_map_to_list()
        
        self.bitString = bitString


        self.build_ansatz()
        
    def build_ansatz(self):

        if self.type == 'two_local':
            return self.two_local()
        
        if self.type == 'heuristic':
            if self.bitString[0]:
                return self.heuristic_bit()
            else:
                return self.heuristic()

    def two_local (self):
        return TwoLocal(num_qubits=self.n_qubits, rotation_blocks=[self.rotation_blocks],
                                          entanglement_blocks=self.entanglement_blocks, entanglement=self.entangler_map, skip_unentangled_qubits=self.reduce_parameters,
                                          reps=self.depth, parameter_prefix='', insert_barriers=True, initial_state=self.initial_state)
    
    def _old_heuristic (self):

        self.var_form = QuantumCircuit(self.n_qubits)
        if self.initial_state == 'HF' and bool(self.starting_occupations):
            if type(self.starting_occupations) is int:
                t_list = [0]*self.n_qubits
                num_orb = int(self.n_qubits/2)
                i = -1
                spin_variable = int(self.starting_occupations//2) #beta
                while i <spin_variable:
                    i+=1
                    t_list[i+num_orb]+=1
                i = -1
                spin_variable = self.starting_occupations - spin_variable #alpha
                while i <spin_variable:
                    i+=1
                    t_list[i]+=1
                self.starting_occupations = t_list[::-1] #in quantum library in little endian (a.k.a. qiskit)
            for idx, i in enumerate(self.starting_occupations):
                if i:
                    self.var_form.x(idx)
            #self.var_form = HartreeFock()
            
        self.string = '-1'

        
        list_ent = []
        for index in range(0,len(self.entangler_map)):
            list_ent.append(self.entangler_map[index])
        rep = 1
        cycle = 1

        #-----ITERATION VARIATIONAL LAYER-----
        while cycle <= self.depth:
            if rep > len(self.entangler_map):
                rep = rep - len(self.entangler_map)
            ent_repetition = list_ent[rep-1]


            #-----ROTATIONS FOR i-th LAYER------
            for r_gate in self.rotation_blocks:
                for i in range(0,self.n_qubits):
                    #TODO: se abbiamo mappe diverse per layer diversi come si comporta la parametrizzazione ridotta?
                    if self.reduce_parameters:
                        if i not in self.entangled_qubit_list:
                            continue
                    self.string = str(int(self.string) + 1)
                    getattr(self.var_form, r_gate)(Parameter(self.string), qubit=i)  

            #-----ENTANGLING FOR i-th LAYER----
            for c,t in ent_repetition:
                if self.entanglement_blocks in _PARAM_ENT_GATE_LIST:
                    self.string = str(int(self.string) + 1)
                    getattr(self.var_form, self.entanglement_blocks)(Parameter(self.string), c, t)

                elif self.entanglement_blocks in _ENT_GATE_LIST:
                    getattr(self.var_form, self.entanglement_blocks)(c, t)
            self.var_form.barrier()

            rep = rep+1    
            cycle = cycle+1
        
        #-----LAST ROTATION BLOCK------------
        for r_gate in self.rotation_blocks:
            for i in range(0,self.n_qubits):
                self.string = str(int(self.string) + 1)
                getattr(self.var_form, r_gate)(Parameter(self.string), qubit=i)

        

            #TODO: qubit singolarmente parametrizzati dovrebbero avere il parametro nell'ultimo layer così da evitare la decoerenza

#------------------------------------------------------------
#--------------------V-SHAPE ADDITIONAL LAYER----------------
#-----------------------(Bitstring version)------------------
#------------------------------------------------------------

    def add_bitLayer_og(self, additional_entangler_map = None, repetition = 0, fixed_parameters=None):
        if additional_entangler_map == None:
            additional_entangler_map = self.additional_entangler_map


        maxCorrelation = 10#int(np.ceil(np.log2((max(self.depth)+1)*(self.n_qubits*len(self.rotation_blocks))))) + 1
        self.counter = len(self.var_form.parameters) - 1
        
        if len(fixed_parameters)>0:
            self.var_form.assign_parameters(fixed_parameters, True)  
            self.var_form.barrier()
            self.counter = len(fixed_parameters) - 1
        list_ent = []
        for index in range(0,len(additional_entangler_map)):
            list_ent.append(additional_entangler_map[index])
        
        
        for i in range(0,len(list_ent)):
            flatten_ = set(np.array(list_ent[i]).flatten())
            for rep in range(repetition):
                #ROTATIONS LAYER
                for s in flatten_:
                    self.counter = self.counter + 1
                    self.var_form.ry(Parameter(parameterBitString(maxCorrelation,self.counter)),s)
                #DESCENDING ENTANGLING LADDER V-shape
                for c,t in list_ent[i]:
                    if self.entanglement_blocks == 'su4':
                        self.general_real_SU4(c,t,maxCorrelation)
                    else:
                        self.var_form.cx(c,t)

                #CENTRAL ROTATIONS LAYER
                self.var_form.barrier()
                for s in flatten_:
                    self.counter = self.counter + 1
                    self.var_form.ry(Parameter(parameterBitString(maxCorrelation,self.counter)),s)
                self.var_form.barrier()
                #ASCENDING ENTANGLING LADDER V-shape
                list_ent_rev = list_ent[i][::-1]
                for c,t in list_ent_rev:
                    if self.entanglement_blocks == 'su4':
                        self.general_real_SU4(c,t,maxCorrelation)
                    else:
                        self.var_form.cx(c,t)
                self.var_form.barrier()

            #FINAL ROTATIONS LAYER
            for s in flatten_:
                self.counter = self.counter + 1
                self.var_form.ry(Parameter(parameterBitString(maxCorrelation,self.counter)),s)
        return() 
    


    def add_bitLayer_test(self, additional_entangler_map = None, repetition = 0, fixed_parameters=None, entangling_gate = None):
        if entangling_gate == None:
            entangling_gate = self.entanglement_blocks
        if additional_entangler_map == None:
            additional_entangler_map = self.additional_entangler_map


        maxCorrelation = 10#int(np.ceil(np.log2((max(self.depth)+1)*(self.n_qubits*len(self.rotation_blocks))))) + 1
        self.counter = len(self.var_form.parameters) - 1
        
        if len(fixed_parameters)>0:
            self.var_form.assign_parameters(fixed_parameters, True)  
            self.var_form.barrier()
            self.counter = len(fixed_parameters) - 1
        list_ent = []
        for index in range(0,len(additional_entangler_map)):
            list_ent.append(additional_entangler_map[index])
        
        
        for i in range(0,len(list_ent)):
            flatten_ = set(np.array(list_ent[i]).flatten())
            for rep in range(repetition):
                #ROTATIONS LAYER
                for s in flatten_:
                    self.counter = self.counter + 1
                    self.var_form.ry(Parameter(parameterBitString(maxCorrelation,self.counter)),s)
                #DESCENDING ENTANGLING LADDER V-shape
                for c,t in list_ent[i]:
                    if entangling_gate in _PARAM_ENT_GATE_LIST:
                        if entangling_gate == 'su4':
                            self.general_real_SU4(c, t, self.maxCorrelation)
                        else:   

                            self.counter = self.counter + 1
                            getattr(self.var_form, entangling_gate)(Parameter(parameterBitString(self.maxCorrelation,self.counter)), c, t)

                    elif entangling_gate in _ENT_GATE_LIST:
                        getattr(self.var_form, entangling_gate)(c, t)

                #CENTRAL ROTATIONS LAYER
                self.var_form.barrier()
                for s in flatten_:
                    self.counter = self.counter + 1
                    self.var_form.ry(Parameter(parameterBitString(maxCorrelation,self.counter)),s)
                self.var_form.barrier()
                #ASCENDING ENTANGLING LADDER V-shape
                list_ent_rev = list_ent[i][::-1]
                for c,t in list_ent_rev:
                    if entangling_gate in _PARAM_ENT_GATE_LIST:
                        if entangling_gate == 'su4':
                            self.general_real_SU4(c, t, self.maxCorrelation)
                        else:   

                            self.counter = self.counter + 1
                            getattr(self.var_form, entangling_gate)(Parameter(parameterBitString(self.maxCorrelation,self.counter)), c, t)

                    elif entangling_gate in _ENT_GATE_LIST:
                        getattr(self.var_form, entangling_gate)(c, t)
                self.var_form.barrier()

            #FINAL ROTATIONS LAYER
            for s in flatten_:
                self.counter = self.counter + 1
                self.var_form.ry(Parameter(parameterBitString(maxCorrelation,self.counter)),s)
        return() 
    
    #TESTING FOR BITSTRIN ENUMERATION
    def heuristic_bit (self):

        self.var_form = QuantumCircuit(self.n_qubits)
        if self.initial_state == 'HF' and bool(self.starting_occupations):
            if type(self.starting_occupations) is int:
                t_list = [0]*self.n_qubits
                num_orb = int(self.n_qubits/2)
                i = -1
                spin_variable = int(self.starting_occupations//2) #beta
                while i <spin_variable:
                    i+=1
                    t_list[i+num_orb]+=1
                i = -1
                spin_variable = self.starting_occupations - spin_variable #alpha
                while i <spin_variable:
                    i+=1
                    t_list[i]+=1
                self.starting_occupations = t_list[::-1] #in quantum library in little endian (a.k.a. qiskit)
            for idx, i in enumerate(self.starting_occupations):
                if i:
                    self.var_form.x(idx)
            #self.var_form = HartreeFock()

        # COMPUTE THE LENGHT OF THE BITSTRING
        self.maxCorrelation = 10#int(np.ceil(np.log2((max(self.depth)+1)*(self.n_qubits*len(self.rotation_blocks))))) + 1
        # IT IS POSSIBLE TO DEFINE A CIRCUIT WITH AN OFFSET APPLIED ON THE ENUMERATION OF THE PARAMETERS
        self.counter = -1 + self.bitString[1] 

        #--------ALTERNATED MAP REPETITION---------
        if len(self.depth)==1 and len(self.entangler_map)>1:
            list_ent = []
            for index in range(0,len(self.entangler_map)):
                list_ent.append(self.entangler_map[index])
            rep = 1
            cycle = 1

            #-----ITERATION VARIATIONAL LAYER-----
            while cycle <= self.depth[0]:
                if rep > len(self.entangler_map):
                    rep = rep - len(self.entangler_map)
                ent_repetition = list_ent[rep-1]


                #-----ROTATIONS FOR i-th LAYER------
                for r_gate in self.rotation_blocks:
                    for i in range(0,self.n_qubits):
                        #TODO: se abbiamo mappe diverse per layer diversi come si comporta la parametrizzazione ridotta?
                        if self.reduce_parameters:
                            if i not in self.entangled_qubit_list:
                                continue
                        self.counter = self.counter + 1
                        getattr(self.var_form, r_gate)(Parameter(parameterBitString(self.maxCorrelation,self.counter)), qubit=i)  

                #-----ENTANGLING FOR i-th LAYER----
                for c,t in ent_repetition:
                    if self.entanglement_blocks in _PARAM_ENT_GATE_LIST:
                        if self.entanglement_blocks == 'su4':
                            self.general_real_SU4(c, t, self.maxCorrelation)
                        else:   

                            self.counter = self.counter + 1
                            getattr(self.var_form, self.entanglement_blocks)(Parameter(parameterBitString(self.maxCorrelation,self.counter)), c, t)

                    elif self.entanglement_blocks in _ENT_GATE_LIST:
                        getattr(self.var_form, self.entanglement_blocks)(c, t)
                self.var_form.barrier()

                rep = rep+1    
                cycle = cycle+1
            
            #-----LAST ROTATION BLOCK------------
            for r_gate in self.rotation_blocks:
                for i in range(0,self.n_qubits):
                    self.counter = self.counter + 1 
                    getattr(self.var_form, r_gate)(Parameter(parameterBitString(self.maxCorrelation,self.counter)), qubit=i)

        
        else :
            rep = 0
            cycle = 0

            for map_depth in self.depth:
                ent_repetition = self.entangler_map[rep]
                for i in range(map_depth):
                    
                    #-----ROTATIONS FOR i-th LAYER------
                    for r_gate in self.rotation_blocks:
                        for i in range(0,self.n_qubits):
                            #TODO: se abbiamo mappe diverse per layer diversi come si comporta la parametrizzazione ridotta?
                            if self.reduce_parameters:
                                if i not in self.entangled_qubit_list:
                                    continue
                            self.counter = self.counter + 1
                            getattr(self.var_form, r_gate)(Parameter(parameterBitString(self.maxCorrelation,self.counter)), qubit=i)  
                                
                    #-----ENTANGLING FOR i-th LAYER----
                    for c,t in ent_repetition:
                        if self.entanglement_blocks in _PARAM_ENT_GATE_LIST:

                            if self.entanglement_blocks == 'su4':
                                self.general_real_SU4(c, t,10)
                            else:
                                self.counter = self.counter + 1
                                getattr(self.var_form, self.entanglement_blocks)(Parameter(parameterBitString(self.maxCorrelation,self.counter)), c, t)

                        elif self.entanglement_blocks in _ENT_GATE_LIST:
                            getattr(self.var_form, self.entanglement_blocks)(c, t)
                    self.var_form.barrier()
                    cycle = cycle+1
                rep = rep+1
                if cycle == np.sum(self.depth):
                    #-----LAST ROTATION BLOCK------------
                    for r_gate in self.rotation_blocks:
                        for i in range(0,self.n_qubits):
                            self.counter = self.counter + 1
                            getattr(self.var_form, r_gate)(Parameter(parameterBitString(self.maxCorrelation,self.counter)), qubit=i)
    


    def general_SU4(self, q0,q1,M):
        def add_single_SU2 (self,q, M):

            self.counter = self.counter + 1
            self.var_form.rz(Parameter(parameterBitString(M,self.counter)),q)

            self.counter = self.counter + 1
            self.var_form.ry(Parameter(parameterBitString(M,self.counter)),q)

            self.counter = self.counter + 1
            self.var_form.rz(Parameter(parameterBitString(M,self.counter)),q)


        def add_N_block(self, q0,q1, M):
            self.var_form.cx(q1, q0)

            self.counter = self.counter + 1
            self.var_form.rz(Parameter(parameterBitString(M,self.counter)),q0)


            self.counter = self.counter + 1
            self.var_form.ry(Parameter(parameterBitString(M,self.counter)),q1)

            self.var_form.cx(q0, q1)


            self.counter = self.counter + 1
            self.var_form.ry(Parameter(parameterBitString(M,self.counter)),q1)

            self.var_form.cx(q1, q0)

        
        add_single_SU2(self, q0, M)
        add_single_SU2(self, q1, M)
        self.var_form.barrier()
        add_N_block(self, q0,q1, M)
        self.var_form.barrier()
        add_single_SU2(self, q0, M)
        add_single_SU2(self, q1, M)
        self.var_form.barrier()


    def general_real_SU4(self, q0,q1,M):
        def add_single_SU2 (self,q, M):

            self.counter = self.counter + 1
            self.var_form.ry(Parameter(parameterBitString(M,self.counter)),q)

        def add_N_block(self, q0,q1, M):
            self.var_form.cx(q1, q0)

            self.counter = self.counter + 1
            self.var_form.ry(Parameter(parameterBitString(M,self.counter)),q1)

            self.var_form.cx(q0, q1)

            self.counter = self.counter + 1
            self.var_form.ry(Parameter(parameterBitString(M,self.counter)),q1)

            self.var_form.cx(q1, q0)

        
        add_single_SU2(self, q0, M)
        add_single_SU2(self, q1, M)
        self.var_form.barrier()
        add_N_block(self, q0,q1, M)
        self.var_form.barrier()
        add_single_SU2(self, q0, M)
        add_single_SU2(self, q1, M)
        self.var_form.barrier()






    def heuristic (self):

        self.var_form = QuantumCircuit(self.n_qubits)
        if self.initial_state == 'HF' and bool(self.starting_occupations):
            if type(self.starting_occupations) is int:
                t_list = [0]*self.n_qubits
                num_orb = int(self.n_qubits/2)
                i = -1
                spin_variable = int(self.starting_occupations//2) #beta
                while i <spin_variable:
                    i+=1
                    t_list[i+num_orb]+=1
                i = -1
                spin_variable = self.starting_occupations - spin_variable #alpha
                while i <spin_variable:
                    i+=1
                    t_list[i]+=1
                self.starting_occupations = t_list[::-1] #in quantum library in little endian (a.k.a. qiskit)
            for idx, i in enumerate(self.starting_occupations):
                if i:
                    self.var_form.x(idx)
            #self.var_form = HartreeFock()
            

        self.string = '-1'

        #--------ALTERNATED MAP REPETITION---------
        if len(self.depth)==1 and len(self.entangler_map)>1:
            list_ent = []
            for index in range(0,len(self.entangler_map)):
                list_ent.append(self.entangler_map[index])
            rep = 1
            cycle = 1

            #-----ITERATION VARIATIONAL LAYER-----
            while cycle <= self.depth[0]:
                if rep > len(self.entangler_map):
                    rep = rep - len(self.entangler_map)
                ent_repetition = list_ent[rep-1]


                #-----ROTATIONS FOR i-th LAYER------
                for r_gate in self.rotation_blocks:
                    for i in range(0,self.n_qubits):
                        #TODO: se abbiamo mappe diverse per layer diversi come si comporta la parametrizzazione ridotta?
                        if self.reduce_parameters:
                            if i not in self.entangled_qubit_list:
                                continue
                        self.string = str(int(self.string) + 1)
                        getattr(self.var_form, r_gate)(Parameter(self.string), qubit=i)  

                #-----ENTANGLING FOR i-th LAYER----
                for c,t in ent_repetition:
                    if self.entanglement_blocks in _PARAM_ENT_GATE_LIST:
                        self.string = str(int(self.string) + 1)
                        getattr(self.var_form, self.entanglement_blocks)(Parameter(self.string), c, t)

                    elif self.entanglement_blocks in _ENT_GATE_LIST:
                        getattr(self.var_form, self.entanglement_blocks)(c, t)
                self.var_form.barrier()

                rep = rep+1    
                cycle = cycle+1
            
            #-----LAST ROTATION BLOCK------------
            for r_gate in self.rotation_blocks:
                for i in range(0,self.n_qubits):
                    self.string = str(int(self.string) + 1)
                    getattr(self.var_form, r_gate)(Parameter(self.string), qubit=i)

        
        else :
            rep = 0
            cycle = 0

            for map_depth in self.depth:
                ent_repetition = self.entangler_map[rep]
                for i in range(map_depth):
                    
                    #-----ROTATIONS FOR i-th LAYER------
                    for r_gate in self.rotation_blocks:
                        for i in range(0,self.n_qubits):
                            #TODO: se abbiamo mappe diverse per layer diversi come si comporta la parametrizzazione ridotta?
                            if self.reduce_parameters:
                                if i not in self.entangled_qubit_list:
                                    continue
                            self.string = str(int(self.string) + 1)
                            getattr(self.var_form, r_gate)(Parameter(self.string), qubit=i)  
                                
                    #-----ENTANGLING FOR i-th LAYER----
                    for c,t in ent_repetition:
                        if self.entanglement_blocks in _PARAM_ENT_GATE_LIST:
                            self.string = str(int(self.string) + 1)
                            getattr(self.var_form, self.entanglement_blocks)(Parameter(self.string), c, t)

                        elif self.entanglement_blocks in _ENT_GATE_LIST:
                            getattr(self.var_form, self.entanglement_blocks)(c, t)
                    self.var_form.barrier()
                    cycle = cycle+1
                rep = rep+1
                if cycle == np.sum(self.depth):
                    #-----LAST ROTATION BLOCK------------
                    for r_gate in self.rotation_blocks:
                        for i in range(0,self.n_qubits):
                            self.string = str(int(self.string) + 1)
                            getattr(self.var_form, r_gate)(Parameter(self.string), qubit=i)

    def _entangler_map_to_list(self):

        """
            Helper function that converts the entangler_map into a simple list of qubits that are involved in the entangling process.
        """ 
        entangled_qubit = []
        for circuit in self.entangler_map:
            for layer in circuit:
                t, c = layer
                if t not in entangled_qubit: entangled_qubit.append(t)
                if c not in entangled_qubit: entangled_qubit.append(c)

        self.entangled_qubit_list = entangled_qubit

def parameterBitString(M,i):
    return format(i,'b').zfill(M)

def depth_given_qubitnumber(n_qubits):
    return n_qubits//2 +1


def separable_ansatz_constructor(input, entangler_maps : list = None, depth : list = None, classical_threshold : int = 4) -> dict:
    internal_input = input
    if entangler_maps is None:
        internal_entangler_maps = internal_input.ansatz_entangler_map
    else:
        internal_entangler_maps = entangler_maps
    
    if depth is None:
        internal_depth = internal_input.ansatz_depth
    else:
        internal_depth = depth
    #if len(np.shape(internal_depth)) == 1:     #an attempt o f generiled depth input
    #    internal_depth = [internal_depth for _ in range(len(internal_entangler_maps))]

    internal_ansatz_origin = len(np.shape(internal_entangler_maps))
    internal_circuits = {}
    if internal_ansatz_origin == 2: #QIDA output
        G = nx.Graph(entangler_maps)
        internal_qubits_subgroups = list(nx.connected_components(G))
        if len(internal_qubits_subgroups) == 1:
            warnings.warn("The graph is fully connected, the ansatz is not separable!")
        for qubits_group in internal_qubits_subgroups:
            local_ansatz = []
            for pair in internal_entangler_maps:
                involved_qubits = set(pair)
                if involved_qubits.issubset(qubits_group):
                    local_ansatz.append([qubits_group.index(pair[0]), qubits_group.index(pair[1])])
            local_num_qubits = len(qubits_group)
            if local_num_qubits > classical_threshold:
                internal_circuits[qubits_group] = {
                    'ansatz' : Ansatz(input, n_qubits=local_num_qubits, ansatz_type = 'heuristic', entangler_map = [local_ansatz], depth = [depth_given_qubitnumber(local_num_qubits)]),
                    'operator' : None,
                    'operator_contribution' : None,
                    'last_results' : None,
                    'last_parameters' : None,
                    'psiopsi' : None,
                    'operator_exp_value' : None,
                    'convergence' : False,
                    'device' : 'QPU'
                    }
            else:
                internal_circuits[qubits_group] = {
                    'operator' : None,
                    'operator_contribution' : None,
                    'last_results' : None,
                    'operator_exp_value' : None,
                    'device' : 'CPU'
                }
    elif internal_ansatz_origin == 4:#input file
        internal_qubits_subgroups = internal_entangler_maps
        if len(internal_qubits_subgroups) == 1:
            warnings.warn("The graph is fully connected, the ansatz is not separable!, are you sure about the input file?")
        for ansatz_idx, ansatz in enumerate(internal_entangler_maps):
            ansatz_key = set()
            for ansatz_piece in ansatz:
                for qubit_pair in ansatz_piece:
                    ansatz_key |= set(qubit_pair)
            local_num_qubits = len(ansatz_key)
            internal_circuits[ansatz_key] = {
                'operator' : None,
                'operator_contribution' : None,
                'last_results' : None,
                'last_parameters' : None,
                'operator_exp_value' : None
                }
            if local_num_qubits > classical_threshold:
                internal_circuits[ansatz_key]['psiopsi'] = None
                internal_circuits[ansatz_key]['ansatz'] = Ansatz(input, n_qubits=local_num_qubits, ansatz_type = 'heuristic', entangler_map = ansatz, depth = internal_depth[ansatz_idx])
                internal_circuits[ansatz_key]['convergence'] = False
                internal_circuits[ansatz_key]['device'] = 'QPU'
            else:
                internal_circuits[ansatz_key]['device'] = 'CPU'
                ## TODO: on Separable_VQE you must contruct a local operator diagonalization
    else:
        raise ValueError("The ansatz is not in a valid format!")

    return internal_circuits