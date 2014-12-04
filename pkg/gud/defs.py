__all__ = ['exfflds', 'checkflds', 'checklflds']

import fortran

# add here new fields to be read in through namelist gud_forcing_params
# first entry is name of the variable to read into
# second is the prefix for runtime paramaters such as PARfile
exfflds = [
    ['surfPAR', 'PAR'],
    ['inputFe', 'iron'],
    ['iceFrac', 'ice'],
    ['windSpeed', 'wind'],
    ['atmospCO2', 'pCO2']
    ]

#conditional = {
#    'alk':  ('#ifdef ALLOW_CARBON\n', '\n#endif'),
#    'o2':   ('#ifdef ALLOW_CARBON\n', '\n#endif'),
#    'cdom': ('#ifdef ALLOW_CDOM\n', '\n#endif'),
#}

p,conds = fortran.readparameters('GUD_SIZE.h', 'GUD_INDICES.h', conditions=True)
conditional = {}
for k,v in conds.items():
    v = [c for c in v if not c.endswith('ALLOW_GUD\n')]
    if k[0] == 'i' and len(v):
        conditional[k[1:]] = (''.join(v), len(v)*'#endif\n')

checkflds = []
checklflds = []
for k in p:
    if k[0] == 'i' and k[1:4] not in ['Min', 'Max']:
        name = k[1:]
        pre,suf = conditional.get(name.lower(), ('', ''))
        if 'e'+name in p:
            checklflds.append((name, '{0:<5s}'.format(name), pre, suf))
        else:
            checkflds.append((name, '{0:<6s}'.format(name), pre, suf))

