class Typical_Bond:
    def __init__(self, length, identity, specific_id):
        self.length = length
        self.identity = identity
        self.specific_id = specific_id
        
#CC_sg = Typical_Bond(1.54, 'CC', 'CC_sg')
#CC_db = Typical_Bond(1.34, 'CC', 'CC_db')


#carbonyl_RC = Typical_Bond(1.52, 'CC', 'X-GR', 'carbonyl')
#carboxyl_RC = Typical_Bond(1.52, 'CC', 'X-GR', 'carboxyl')
#Aro_C = Typical_Bond(1.52, 'CC', 'X-GR', 'aroc')

###LATER### epoxy_RC = Typical_Bond(1.52, 'CC', 'X-GR', 'epoxy')
###LATER### epoxy_CC = Typical_Bond(1.54, 'CC', 'GR-GR', 'epoxy')



#CH_simple = Typical_Bond(1.09, 'CH', 'GR-GR', 'chsimple')
#GRA_C = Typical_Bond(1.42, 'CC', 'X-X', 'graphene')

###later### carbonyl_CO = Typical_Bond(1.23, 'CO', 'GR-GR')

#CO_db = Typical_Bond(1.21, 'CO', 'GR-GR')
#carboxyl_CO_db = Typical_Bond(1.20, 'CO', 'GR-GR')
#carboxyl_CO_sg = Typical_Bond(1.34, 'CO', 'GR-GR')

#carboxyl_OH = Typical_Bond(0.98, 'OH', 'GR-GR')

#hydroxyl_CO = Typical_Bond(1.42, 'CO', 'X-GR')
#hydroxyl_OH = Typical_Bond(0.97, 'OH', 'GR-GR')

###Later### epoxy_CO = Typical_Bond(1.42, 'CO', 'GR-GR')
###later#### carbonyl_RC = Typical_Bond(1.52, 'CC', 'X-GR')

#bond_list = [ CC_db,  CH_simple,  carboxyl_CO_sg,  carboxyl_CO_db,  hydroxyl_CO, hydroxyl_OH, epoxy_RC, epoxy_CO,  GRA_C]
# Aro_C,carboxyl_RC,carbonyl_RC,epoxy_CC, CC_sg, carboxyl_OH, CO_db, , carbonyl_CO,

# https://arxiv.org/pdf/1102.3797.pdf
# epoxy_CO 1.46 and the res. C-C 1.51
        

#Carbon-Carbon
#CC_sg = Typical_Bond(1.54, 'CC', 'GR-GR')
CC_db = Typical_Bond(1.34, 'CC', 'GR-GR')
Aro_C = Typical_Bond(1.52, 'CC', 'X-GR') # + CC_sg, + carboxyl_RC
#carboxyl_RC = Typical_Bond(1.52, 'CC', 'X-GR')
GRA_C = Typical_Bond(1.42, 'CC', 'X-X')

#Carbon-Oxygen
#carboxyl_CO_db = Typical_Bond(1.20, 'CO', 'GR-GR')
CO_db = Typical_Bond(1.21, 'CO', 'GR-GR') #carboxyl_CO_db
carboxyl_CO_sg = Typical_Bond(1.34, 'CO', 'GR-GR')
#hydroxyl_CO = Typical_Bond(1.49, 'CO', 'X-GR') #was 1.42
epoxy_CO = Typical_Bond(1.46, 'CO', 'X-GR')

#Carbon-Hydrogen
CH_simple = Typical_Bond(1.09, 'CH', 'GR-GR', )

#Oxygen-Hydrogen
carboxyl_OH = Typical_Bond(0.98, 'OH', 'GR-GR') #same as hydroxyl_OH
#hydroxyl_OH = Typical_Bond(0.97, 'OH', 'GR-GR')

#Nitrogen-Carbon
nitro_C = Typical_Bond(1.49, 'CN', 'X-GR')

bond_list = [CC_db, Aro_C, GRA_C, CO_db, carboxyl_CO_sg, epoxy_CO, CH_simple, carboxyl_OH, nitro_C] # hydroxyl_CO
