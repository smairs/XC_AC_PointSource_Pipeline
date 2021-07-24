import numpy as np
import pickle


def get_vars(region,wave):
    if wave == '450':
        if region == 'IC348':
            known_vars = [1]
            var_names  = ['IC 348 1']
        if region == 'NGC1333':
            known_vars = [0,10,1,4]
            var_names  = ['IRAS 4A','BOLO 40','NoName','NoName']
        if region == 'NGC2024':
            known_vars = []
            var_names  = []
        if region == 'NGC2071':
            known_vars = [1,2]
            var_names  = ['HOPS 358','HOPS 373']
        if region == 'OMC23':
            known_vars = [1,12,29,100]
            var_names  = ['HOPS 88','HOPS 370','HOPS 383','V1017 ORI']
        if region == 'OPHCORE':
            known_vars = [9]
            var_names  = ['OPH 162624']
        if region == 'SERPM':
            known_vars = [0,10,11]
            var_names  = ['SMM 1','EC 53','SMM 10']
        if region == 'SERPS':
            known_vars = [12]
            var_names  = ['IRAS 18270']
        if region == 'DR21C':
            known_vars = []
            var_names = []
        if region == 'DR21N':
            known_vars = []
            var_names = []
        if region == 'DR21S':
            known_vars = []
            var_names = []
        if region == 'M17':
            known_vars = []
            var_names = []
        if region == 'M17SWex':
            known_vars = []
            var_names = []
        if region == 'S255':
            known_vars = []
            var_names = []

    if wave == '850':
        if region == 'IC348':
            known_vars = [1]
            var_names  = ['IC 348 1']
        if region == 'NGC1333':
            known_vars = [0,19,24,1,2]
            var_names  = ['IRAS 4A','BOLO 40','NoName','NoName','NoName']
        if region == 'NGC2024':
            known_vars = []
            var_names  = []
        if region == 'NGC2071':
            known_vars = [1,12]
            var_names  = ['HOPS 358','HOPS 373']
        if region == 'OMC23':
            known_vars = [40,1,57,97]
            var_names  = ['HOPS 88','HOPS 370','HOPS 383','V1017 ORI']
        if region == 'OPHCORE':
            known_vars = [28]
            var_names  = ['OPH 162624']
        if region == 'SERPM':
            known_vars = [0,2,3]
            var_names  = ['SMM 1','EC 53','SMM 10']
        if region == 'SERPS':
            known_vars = [46]
            var_names  = ['IRAS 18270']
        if region == 'DR21C':
            known_vars = []
            var_names = []
        if region == 'DR21N':
            known_vars = []
            var_names = []
        if region == 'DR21S':
            known_vars = []
            var_names = []
        if region == 'M17':
            known_vars = []
            var_names = []
        if region == 'M17SWex':
            known_vars = []
            var_names = []
        if region == 'S255':
            known_vars = []
            var_names = []
    return(known_vars,var_names)
