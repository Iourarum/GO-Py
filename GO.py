#C 10 O 1 (OH) 1 (COOH) 0.5
#C–C bond lengths of 1.42 Å
import numpy as np
import math
from operator import itemgetter, attrgetter
from retypicalbonds import *
import random
import string
from scipy import spatial
import time

def timefunction(f):
    def f_timer(*args, **kwargs):
        with open(timefile1, 'a') as timefile:
            start = time.time()
            result = f(*args, **kwargs)
            end = time.time()
            if ( str(f.__name__) in func_timedict.keys()):
                func_timedict[str(f.__name__)] += float(end-start)
            else:
                func_timedict[str(f.__name__)] = float(end-start)
            if ( str(f.__name__) in func_calldict.keys()):
                func_calldict[str(f.__name__)] += 1
            else:
                func_calldict[str(f.__name__)] = 1
            timefile.write( str(f.__name__) + " " + str(float(end-start)) + " s" + "\n" )
        return result
    return f_timer

class Atom:
    def __init__(self, position, atom_name, residue_name, residue_type_no, x, y, z):
        self.position = position
        self.atom_name = atom_name
        self.residue_name = residue_name
        self.residue_type_no = residue_type_no
        self.x = x
        self.y = y
        self.z = z
        
    @timefunction        
    def links(self, atom_list):
        d, c = distance_list_sort_and_trimm(build_distance_list(self, atom_list))
        if (str(c) =='[2]'):
            return [(self, 0)]
        else: 
            #return distance_list_sort_and_trimm(build_distance_list(self, atom_list))[0]
            return d
        
    @timefunction
    def linked_to_functional_group(self, atom_list):
        return distance_list_sort_and_trimm(build_distance_list(self, atom_list))[1]

    @timefunction
    def positioning(self, atom_list): # should return either inner or outer
        if (len(self.links(atoms)) <= 2): 
            if (str(self.linked_to_functional_group(atoms)) == '[0]'):
                return "edgy"
            elif ( (str(self.linked_to_functional_group(atoms)) == '[1]' )):
                return "busy"
        elif (len(self.links(atoms)) == 3 ): 
            if (str(self.linked_to_functional_group(atoms)) == '[0]'):
                return "center"
    
    @timefunction
    def is_graphene(self):
        if (self.residue_name=='GGG'):
            return "X" 
        else:
            return "GR"    

def show(func_timedict, func_calldict):
    for key in func_timedict:
        print(str(key) + ": " + str(func_timedict[key]) + " s   " + str(func_calldict[key]) + " calls.")
    

@timefunction   
def provide_identity(initial_len, atom_list):
    identity_string = ''
    for i in range(initial_len,len(atom_list)): #iterate the elements in List_B, look up the values in the dictionary, and append them to the list. Also, don't forget to update the dictionary with the new (higher) percentages.
        identity_string = identity_string + str(atom_list[i].atom_name)[0:2] + str(atom_list[i].x)[0:4] + str(atom_list[i].y)[0:4] + str(atom_list[i].z)[0:4]+" "
    return identity_string

@timefunction
def check_identity(given_identity):
    new = 1
    with open(identity_file, 'r+') as f:
        content = f.readlines()
        identities = [x.split() for x in content]
        for identity in identities:
            score = 0
            suma = 0
            if ( len(identity) == len(given_identity)):
                for element in identity:
                    for i in given_identity:
                        if (identities[identity][element] == given_identity[i]):
                            suma += 1
                score = suma / len(given_identity)
            if (score >= 0.95):
                new = 0
        if (new == 1):
            print("should write to file")
            f.write(given_identity + " \n")
    return new


@timefunction
def read_in_graphene(filename1):
    with open(filename1) as f:
        content = f.readlines()
        atom_lines = [x.split() for x in content if 'ATOM' in str(x)]  
     #   print(atom_lines)
        atoms = [Atom(int(str(atom_lines[x][1])), str(atom_lines[x][2]), str(atom_lines[x][3]), int(str(atom_lines[x][4])), float(str(atom_lines[x][5])), float(str(atom_lines[x][6])), float(str(atom_lines[x][7]))) for x in range(len(atom_lines))] 
    return atoms

#def using_kdtree2(points1, points2, cutoff):
#    # build the KDTree using the *larger* points array
#    tree = spatial.cKDTree(points1)
#    groups = tree.query_ball_point(points2, cutoff)
#    indices = np.unique(IT.chain.from_iterable(groups))
#    return indices

@timefunction
def calculate_3Ddistance_2atoms(atom1, atom2):
    return spatial.distance.euclidean((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
#    return math.sqrt((atom1.x - atom2.x)**2 + (atom1.y - atom2.y)**2 + (atom1.z - atom2.z)**2)

@timefunction
def build_distance_list(given_atom, atom_list):
    distance_list = [(atom, calculate_3Ddistance_2atoms(given_atom, atom)) for atom in atom_list if ((abs(given_atom.x - atom.x) <= 2) and (abs(given_atom.y - atom.y) <= 2) and (abs(given_atom.z - atom.z) <= 2))]
    return distance_list

@timefunction
def distance_list_sort_and_trimm(distance_list):
    sorted_list = sorted(distance_list, key = itemgetter(1))
    specific_list = [x for x in bond_list if (sorted_list[0][0].atom_name[:1] in x.identity)]
    trimmed_list = [x for x in sorted_list if ( 2 >= x[1])]#(bond_length_max_threshold * (max(specific_list, key = attrgetter('length'))).length >= x[1])]
#    print(trimmed_list[0][0].atom_name)
    del trimmed_list[0]
    identified_bonds = []
    for element in trimmed_list:
        atom_couple_name = sorted(sorted_list[0][0].atom_name[:1] + element[0].atom_name[:1])
        atom_couple_chain = sorted( str(sorted_list[0][0].is_graphene()) + str(element[0].is_graphene())) 
        functional_group_link = [1 if ( 'GR' in str(sorted_list[0][0].is_graphene()) + str(element[0].is_graphene())) else 0]
        for bond in specific_list:
            if (atom_couple_name == sorted(bond.identity)):
                if (bond_length_min_threshold * bond.length <= element[1] <= bond_length_max_threshold * bond.length):
                    identified_bonds.append((element[0], element[1], atom_couple_name, atom_couple_chain)) 
  #  print ('problem: ', len(trimmed_list))
    if (len(trimmed_list) == len(identified_bonds)):
        return identified_bonds, functional_group_link
    else:
        return identified_bonds, [2]

@timefunction
def top_or_down():
    ct = random.randint(1,2)
    if (ct == 1):
        return 1
    else:
        return -1

@timefunction
def availability_scanner_epo(atom_list, neighbours_list):
    scanner = []
    for atom in neighbours_list:
        if (int(str(atom.linked_to_functional_group(atom_list))[1]) == 0):
            scanner.append(atom)
    return scanner

@timefunction
def availability_scanner(atom_list, fgroup):
    scanner = []
    if (fgroup == "epoxy" or fgroup == "hydroxyl"):
        for atom in atom_list:
            if (int(str(atom.linked_to_functional_group(atom_list))[1]) == 0):
                scanner.append(atom)        
    elif (fgroup == "cooh"):
        for atom in atom_list:
            if (atom.positioning(atoms) == "edgy"):
                scanner.append(atom)
    return scanner

@timefunction
def add_carboxyl(atom_pos, x, y, z, ct, atom_list, added_functional_groups, added_functional_groups_COOH, aatom):
    current_size = len(atom_list)
    placed = 0
    alpha = random.randint(0,359)
    while (placed <=359):
        alpha += 1
        carbon_atom = Atom(current_size + 1, 'C4', 'C1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(x)), float("{0:.3f}".format(y)), float("{0:.3f}".format(ct * 1.52 + z)))
        atom_list.append(carbon_atom)
        if ( (carbon_atom.links(atom_list)[0][0].position == atom_pos) and (len(carbon_atom.links(atom_list)) == 1) ):
            h = math.sin(math.radians(60)) * 1.20
            oxygen_atom_1 = Atom(current_size + 2, 'OJ', 'C1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(carbon_atom.x - math.cos(math.radians(alpha)) * h)), float("{0:.3f}".format(carbon_atom.y - math.sin(math.radians(alpha)) * h)), float("{0:.3f}".format(carbon_atom.z + ct * math.cos(math.radians(60)) * 1.20)))
            atom_list.append(oxygen_atom_1)
            if ( (oxygen_atom_1.links(atom_list)[0][0].position == carbon_atom.position) and (len(oxygen_atom_1.links(atom_list)) == 1) ):
                h = math.sin(math.radians(60)) * 1.34
                oxygen_atom_2 = Atom(current_size + 3, 'OK', 'C1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(carbon_atom.x - math.cos(math.radians(alpha + 180)) * h)), float("{0:.3f}".format(carbon_atom.y - math.sin(math.radians(alpha+180)) * h)), float("{0:.3f}".format(carbon_atom.z + ct * math.cos(math.radians(60)) * 1.34)) )
                atom_list.append(oxygen_atom_2)
                if ( (oxygen_atom_2.links(atom_list)[0][0].position == carbon_atom.position) and (len(oxygen_atom_2.links(atom_list)) == 1) ):
                    hydrogen_atom = Atom(current_size + 4, 'HJ', 'C1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(oxygen_atom_2.x)),  float("{0:.3f}".format(oxygen_atom_2.y)),  float("{0:.3f}".format(oxygen_atom_2.z+ct*0.98)))
                    atom_list.append(hydrogen_atom)
                    if ( (hydrogen_atom.links(atom_list)[0][0].position == oxygen_atom_2.position) and (len(hydrogen_atom.links(atom_list)) == 1) ):
                        placed = 888
                        link_map_cooh.remove(aatom)
                        if (aatom in link_map_central):
                            link_map_central.remove(aatom)
                    else: 
                        placed += 1
                        del hydrogen_atom
                        del atom_list[current_size + 3]
                        del oxygen_atom_2
                        del atom_list[current_size + 2]
                        del oxygen_atom_1
                        del atom_list[current_size + 1]
                        del carbon_atom
                        del atom_list[current_size + 0]
                else:
                    placed += 1
                    del oxygen_atom_2
                    del atom_list[current_size + 2]
                    del oxygen_atom_1
                    del atom_list[current_size + 1]
                    del carbon_atom
                    del atom_list[current_size + 0]
            else: 
                placed += 1
                del oxygen_atom_1
                del atom_list[current_size + 1]
                del carbon_atom
                del atom_list[current_size + 0]
        else:
            del carbon_atom 
            del atom_list[current_size + 0]
            carbon_atom = Atom(current_size + 1, 'C4', 'C1A', str(added_functional_groups + 1 + 1596), x, y, -1 * ct * 1.52 + z)
            atom_list.append(carbon_atom)
            if ( (carbon_atom.links(atom_list)[0][0].position == atom_pos) and (len(carbon_atom.links(atom_list)) == 1) ):
                h = math.sin(math.radians(60)) * 1.20
                oxygen_atom_1 = Atom(current_size + 2, 'OJ', 'C1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(carbon_atom.x - math.cos(math.radians(alpha)) * h)), float("{0:.3f}".format(carbon_atom.y - math.sin(math.radians(alpha)) * h)), float("{0:.3f}".format(carbon_atom.z + -1 * ct * math.cos(math.radians(60)) * 1.20)))
                atom_list.append(oxygen_atom_1)
                if ( (oxygen_atom_1.links(atom_list)[0][0].position == carbon_atom.position) and (len(oxygen_atom_1.links(atom_list)) == 1) ):
                    h = math.sin(math.radians(60)) * 1.34
                    oxygen_atom_2 = Atom(current_size + 3, 'OK', 'C1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(carbon_atom.x - math.cos(math.radians(alpha + 180)) * h)), float("{0:.3f}".format(carbon_atom.y - math.sin(math.radians(alpha+180)) * h)), float("{0:.3f}".format(carbon_atom.z + -1 * ct * math.cos(math.radians(60)) * 1.34)) )
                    atom_list.append(oxygen_atom_2)
                    if ( (oxygen_atom_2.links(atom_list)[0][0].position == carbon_atom.position) and (len(oxygen_atom_2.links(atom_list)) == 1) ):
                        hydrogen_atom = Atom(current_size + 4, 'HJ', 'C1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(oxygen_atom_2.x)),  float("{0:.3f}".format(oxygen_atom_2.y)),  float("{0:.3f}".format(oxygen_atom_2.z + -1 * ct * 0.98)))
                        atom_list.append(hydrogen_atom)
                        if ( (hydrogen_atom.links(atom_list)[0][0].position == oxygen_atom_2.position) and (len(hydrogen_atom.links(atom_list)) == 1) ):
                            placed = 888
                            link_map_cooh.remove(aatom)
                            if (aatom in link_map_central):
                                link_map_central.remove(aatom)
                        else: 
                            placed += 1
                            del hydrogen_atom
                            del atom_list[current_size + 3]
                            del oxygen_atom_2
                            del atom_list[current_size + 2]
                            del oxygen_atom_1
                            del atom_list[current_size + 1]
                            del carbon_atom
                            del atom_list[current_size + 0]
                    else:
                        placed += 1
                        del oxygen_atom_2
                        del atom_list[current_size + 2]
                        del oxygen_atom_1
                        del atom_list[current_size + 1]
                        del carbon_atom
                        del atom_list[current_size + 0]
                else: 
                    placed += 1
                    del oxygen_atom_1
                    del atom_list[current_size + 1]
                    del carbon_atom
                    del atom_list[current_size + 0]
            else:
                placed += 1
                del carbon_atom 
                del atom_list[current_size + 0]
    return atom_list

@timefunction
def add_epoxy(chosen_atom, ct, atom_list, added_functional_groups, added_functional_groups_epoxy):
    current_size = len(atom_list)
    placed = 0
    neighbours = chosen_atom.links(atom_list)
    list_of_n = [x[0] for x in neighbours if ((x[0].atom_name == 'CX') and (int(str(x[0].linked_to_functional_group(atom_list))[1]) == 0) ) ]
    my_scan = availability_scanner_epo(atom_list, list_of_n)
    if (len(my_scan) != 0):
        atom2 = random.choice(my_scan)
        epoxy_atom = Atom(current_size + 1, 'OE', 'E1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(abs(chosen_atom.x - atom2.x) / 2 + min(chosen_atom.x, atom2.x))), float("{0:.3f}".format(abs(chosen_atom.y - atom2.y) / 2  + min(chosen_atom.y, atom2.y))), float("{0:.3f}".format(ct * 1.46 * math.sin(math.radians(60)))))
        atom_list.append(epoxy_atom)
        if (len(epoxy_atom.links(atom_list)) == 2):
            if ( ((epoxy_atom.links(atom_list)[0][0].position == chosen_atom.position) or (epoxy_atom.links(atom_list)[1][0].position == chosen_atom.position)) and ((epoxy_atom.links(atom_list)[1][0].position == atom2.position) or (epoxy_atom.links(atom_list)[0][0].position == atom2.position))):
                placed = 888
                link_map_central.remove(chosen_atom)
                link_map_central.remove(atom2) #added
                if (chosen_atom in link_map_cooh):
                    link_map_cooh.remove(chosen_atom)
                if (atom2 in link_map_cooh):
                    link_map_cooh.remove(atom2)
            else:
                del epoxy_atom
                del atom_list[current_size + 0]
        else:
            del epoxy_atom
            del atom_list[current_size + 0]
            epoxy_atom = Atom(current_size + 1, 'OE', 'E1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(abs(chosen_atom.x - atom2.x) / 2 + min(chosen_atom.x, atom2.x))), float("{0:.3f}".format(abs(chosen_atom.y - atom2.y) / 2  + min(chosen_atom.y, atom2.y))), float("{0:.3f}".format(-1 * ct * 1.46 * math.sin(math.radians(60)))))
            atom_list.append(epoxy_atom)
            if (len(epoxy_atom.links(atom_list)) == 2):
                if ( ((epoxy_atom.links(atom_list)[0][0].position == chosen_atom.position) or (epoxy_atom.links(atom_list)[1][0].position == chosen_atom.position)) and ((epoxy_atom.links(atom_list)[1][0].position == atom2.position) or (epoxy_atom.links(atom_list)[0][0].position == atom2.position))):
                    placed = 888
                    link_map_central.remove(chosen_atom)
                    if (chosen_atom in link_map_cooh):
                        link_map_cooh.remove(chosen_atom)
                else:
                    del epoxy_atom
                    del atom_list[current_size + 0]
            else:
                del epoxy_atom
                del atom_list[current_size + 0]
    return atom_list
 
@timefunction
def add_hydroxyl(chosen_atom, ct, atom_list, added_functional_groups, added_functional_groups_OH): 
    current_size = len(atom_list)
    placed = 0
    alpha = random.randint(0,359)
    while (placed <= 359):
        alpha += 1
        oxygen_atom = Atom(current_size + 1, 'OL', 'H1A', str(added_functional_groups + 1 + 1596), chosen_atom.x, chosen_atom.y, ct * 1.49 + chosen_atom.z)
        atom_list.append(oxygen_atom)
        if ( (oxygen_atom.links(atom_list)[0][0].position == chosen_atom.position) and (len(oxygen_atom.links(atom_list)) == 1)):
            h = math.sin(math.radians(19)) * 0.98
            h_sp = math.cos(math.radians(19)) * 0.98
            hydrogen_atom = Atom(current_size + 2, 'HK', 'H1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(oxygen_atom.x - math.cos(math.radians(alpha)) * h_sp)), float("{0:.3f}".format(oxygen_atom.y - math.sin(math.radians(alpha)) * h_sp)), float("{0:.3f}".format(oxygen_atom.z + ct * h)))
            atom_list.append(hydrogen_atom)
            if ((hydrogen_atom.links(atom_list)[0][0].position == oxygen_atom.position) and (len (hydrogen_atom.links(atom_list)) == 1) ):
                placed = 888
                link_map_central.remove(chosen_atom)
                if (chosen_atom in link_map_cooh):
                    link_map_cooh.remove(chosen_atom)
            else:
                placed += 1
                del hydrogen_atom
                del atom_list[current_size + 1]
                del oxygen_atom
                del atom_list[current_size + 0]
        else:
            del oxygen_atom
            del atom_list[current_size + 0]
            oxygen_atom = Atom(current_size + 1, 'OL', 'H1A', str(added_functional_groups + 1 + 1596), chosen_atom.x, chosen_atom.y,-1 * ct * 1.49 + chosen_atom.z)
            atom_list.append(oxygen_atom)
            if ( (oxygen_atom.links(atom_list)[0][0].position == chosen_atom.position) and (len(oxygen_atom.links(atom_list)) == 1)):
                h = math.sin(math.radians(19)) * 0.98
                h_sp = math.cos(math.radians(19)) * 0.98
                hydrogen_atom = Atom(current_size + 2, 'HK', 'H1A', str(added_functional_groups + 1 + 1596), float("{0:.3f}".format(oxygen_atom.x - math.cos(math.radians(alpha)) * h_sp)), float("{0:.3f}".format(oxygen_atom.y - math.sin(math.radians(alpha)) * h_sp)), float("{0:.3f}".format(oxygen_atom.z + -1 * ct * h)))
                atom_list.append(hydrogen_atom)
                if ((hydrogen_atom.links(atom_list)[0][0].position == oxygen_atom.position) and (len (hydrogen_atom.links(atom_list)) == 1) ):
                    placed = 888
                    link_map_central.remove(chosen_atom)
                    if (chosen_atom in link_map_cooh):
                        link_map_cooh.remove(chosen_atom)
                else:
                    placed += 1
                    del hydrogen_atom
                    del atom_list[current_size + 1]
                    del oxygen_atom
                    del atom_list[current_size + 0]
            else:
                placed += 1
                del oxygen_atom
                del atom_list[current_size + 0]
    return atom_list

@timefunction
def random_addition_one_carboxyl(atom_list, added_functional_groups, added_functional_groups_COOH):
    temp_list = link_map_cooh
    chosen_atom = random.choice(temp_list)
     #chosen_atom = random.choice(availability_scanner(atom_list, 'cooh'))
    new_atom_list = add_carboxyl(chosen_atom.position, chosen_atom.x, chosen_atom.y, chosen_atom.z, top_or_down(), atom_list, added_functional_groups, added_functional_groups_COOH, chosen_atom)
    return new_atom_list

@timefunction
def random_addition_one_epoxy(atom_list, added_functional_groups, added_functional_groups_epoxy):
    chosen_atom = random.choice(link_map_central)
    #chosen_atom = random.choice(availability_scanner(atom_list, 'epoxy'))
    new_atom_list = add_epoxy(chosen_atom, top_or_down(), atom_list,added_functional_groups, added_functional_groups_epoxy)
    return new_atom_list
    
@timefunction
def random_addition_one_hydroxyl(atom_list, added_functional_groups, added_functional_groups_OH):
    chosen_atom = random.choice(link_map_central)
    #chosen_atom = random.choice(availability_scanner(atom_list, 'hydroxyl'))
    new_atom_list = add_hydroxyl(chosen_atom, top_or_down(), atom_list, added_functional_groups, added_functional_groups_OH)
    return new_atom_list

@timefunction
def create_morphology_COOH(atom_list, no_cooh, no_epoxy, no_hydroxyl, added_functional_groups, initial_len, added_functional_groups_COOH, added_functional_groups_epoxy, added_functional_groups_OH):
    must_add = no_cooh + no_epoxy + no_hydroxyl
    while (must_add > 0):
        print("progress: ", must_add)
        print("no_cooh: ", no_cooh)
        print("no_epoxy: ", no_epoxy)
        print("no_hydroxyl: ", no_hydroxyl)
#        if ((no_cooh != 0) and (no_epoxy != 0) and (no_hydroxyl != 0)):
#            make_list = ['cooh', 'epoxy', 'hydroxyl']
#        elif ((no_cooh == 0) and (no_epoxy != 0) and (no_hydroxyl != 0)):
#            make_list = ['epoxy', 'hydroxyl']
#        elif ((no_cooh == 0) and (no_epoxy == 0) and (no_hydroxyl != 0)):
#            make_list = ['hydroxyl']
#        elif ((no_cooh != 0) and (no_epoxy == 0) and (no_hydroxyl != 0)):
#            make_list = ['cooh', 'hydroxyl']
#        elif ((no_cooh != 0) and (no_epoxy == 0) and (no_hydroxyl == 0)):
#            make_list = ['cooh']
#        elif ((no_cooh != 0) and (no_epoxy != 0) and (no_hydroxyl == 0)):
#            make_list = ['cooh', 'epoxy']
#        elif ((no_cooh == 0) and (no_epoxy != 0) and (no_hydroxyl == 0)):
#            make_list = ['epoxy']
        if (no_cooh != 0):
            make_list = ['cooh']
        elif ((no_epoxy != 0) and (no_hydroxyl != 0) ):
            make_list = ['epoxy', 'hydroxyl']
        elif (no_epoxy != 0):
            make_list = ['epoxy']
        elif (no_hydroxyl != 0):
            make_list = ['hydroxyl']
        chosen = random.choice(make_list)
        print("chosen: ", chosen)
        if (chosen == 'cooh'):
            attempt = 0
            while (attempt < 50):
                old_length = len(atom_list)
                new_atom_list = random_addition_one_carboxyl(atom_list, added_functional_groups, added_functional_groups_COOH)
                if (old_length != len(atom_list)):
                    atom_list = new_atom_list
                    added_functional_groups += 1
                    added_functional_groups_COOH += 1
                    must_add -= 1
                    no_cooh -= 1
                    attempt = 1888
                else:
                    attempt += 1
                print("attempt: ", attempt)
            if (attempt == 50):
                must_add = -1
        elif (chosen == 'epoxy'):
            attempt = 0
            while (attempt < 50):
                old_length = len(atom_list)
                new_atom_list = random_addition_one_epoxy(atom_list, added_functional_groups, added_functional_groups_epoxy)
                if (old_length != len(atom_list)):
                    atom_list = new_atom_list
                    added_functional_groups += 1
                    added_functional_groups_epoxy += 1
                    must_add -= 1
                    no_epoxy -= 1
                    attempt = 1888
                else:
                    attempt += 1
                print("attempt: ", attempt)
            if (attempt == 50):
                must_add = -1
        elif (chosen == 'hydroxyl'):
            attempt = 0
            while (attempt < 50):
                old_length = len(atom_list)
                new_atom_list = random_addition_one_hydroxyl(atom_list, added_functional_groups, added_functional_groups_OH)
                if (old_length != len(atom_list)):
                    atom_list = new_atom_list
                    added_functional_groups += 1
                    added_functional_groups_OH += 1
                    must_add -= 1
                    no_hydroxyl -=1
                    attempt = 1888
                    
                else:
                    attempt += 1
                print("attempt: ", attempt)
            if (attempt == 50):
                must_add = -1
    if (must_add == 0):
        if (check_identity(provide_identity(initial_len, atom_list)) == 1):
            writepdb(atom_list)
        return atom_list, 1
    elif (must_add == -1):
        print("ma = ", -1)
        return atom_list, 0
             
@timefunction       
def lw(max_no, str_obj): #left_whitespaces
    x = max_no - len(str_obj)
    y = 0
    string = ''
    for y in range(x):
        string = string + ' '
    return string

@timefunction
def random_generator(size=8, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

@timefunction
def writepdb(list_of_atoms):
    new = random_generator()
    with open('attempt' + "-" + new + ".pdb", 'a') as le_file:
        for atom in list_of_atoms:
            line = "ATOM" + lw(7, str(atom.position)) + str(atom.position) + lw(4, str(atom.atom_name)) + str(atom.atom_name) + "  " + str(atom.residue_name) + lw(6, str(atom.residue_type_no)) + str(atom.residue_type_no) + lw(12, str(float("{0:.3f}".format(atom.x)))) + str(float("{0:.3f}".format(atom.x))) + lw(8, str(float("{0:.3f}".format(atom.y)))) + str(float("{0:.3f}".format(atom.y))) + lw(8, str(float("{0:.3f}".format(atom.z)))) + str(float("{0:.3f}".format(atom.z))) + "  1.00  0.00             "
            le_file.write(line + '\n')
            print(line)

func_timedict = {}
func_calldict = {}
bond_length_max_threshold = 102.5 / 100
bond_length_min_threshold = 97.5 / 100

#filename1 = r'/home/smuraru/CollectiveData/ScriptsforJuly/JulyScript/RunFolder/new.pdb'
#filename1 = r'/home/smuraru/CollectiveData/ScriptsforJuly/JulyScript/RunFolder/again.pdb'
#filename1 = r'/home/smuraru/CollectiveData/ScriptsforJuly/JulyScript/RunFolder/pristine_graphene.pdb'
#filename1 = r'/home/smuraru/CollectiveData/ScriptsforJuly/JulyScript/RunFolder/newgraphene.pdb'
#filename1 = r'/home/smuraru/ForDecPaper/PG.pdb'
#filename1 = r'/home/smuraru/GOtIt/zigzag-11-6.pdb'
filename1 = r'/home/smuraru/ArtDec/PG_11_5.pdb'
#filename1=r'/home/smuraru/CollectiveData/ScriptsforJuly/JulyScript/RunFolder/1.pdb'
identity_file = r'/home/smuraru/CollectiveData/ScriptsforJuly/JulyScript/RunFolder/identity_file.txt'

timefile1 = r'/home/smuraru/CollectiveData/ScriptsforJuly/JulyScript/RunFolder/timefile.txt'

atoms = read_in_graphene(filename1)

initial_len = len(atoms)
added_functional_groups = 0
added_functional_groups_COOH = 0
added_functional_groups_epoxy = 0
added_functional_groups_OH = 0

link_map_central = availability_scanner(atoms, 'epoxy')
link_map_cooh = availability_scanner(atoms, 'cooh')

new_new_new_atom_list, successful = create_morphology_COOH(atoms, 108, 216, 216, added_functional_groups, initial_len, added_functional_groups_COOH, added_functional_groups_epoxy, added_functional_groups_OH)
