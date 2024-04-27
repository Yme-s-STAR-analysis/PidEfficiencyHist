r'''
    This script is to copy parameter configuration files into StRoot
     > Comment items in [targetes] that you don't need, and run this script
'''

import os

path = '/path/to/params/folder'

targets = {
    'cent': ['CentCorrTool', 'CentParams.h'],
    'dca': ['MeanDcaTool', 'MeanDcaParams.h'],
    'shift': ['TpcShiftTool', 'RunNumber.h'],
    'trigger': ['TriggerTool', 'TriggerTool.cxx'],
    'vtx': ['VtxShiftTool', 'VtxShiftTool.cxx']
}

for k in targets.keys():
    item = targets[k]
    print(f'Now copy item {k}:')
    print(f'-> cp {path}/{item[1]} StRoot/{item[0]}/')
    os.system(f'cp {path}/{item[1]} StRoot/{item[0]}/')

print('All done!')
