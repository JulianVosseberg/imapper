#!/usr/bin/env python

import sys

# Supergroup variables
supergroups6 = {'ACAS' : 'Amoebozoa',
                'EHIS' : 'Amoebozoa',
                'EINV' : 'Amoebozoa',
                'DDIS' : 'Amoebozoa',
                'PPAL' : 'Amoebozoa',
                'ASUB' : 'Amoebozoa',
                'TTRA' : 'Obazoa',
                'SARC' : 'Obazoa',
                'CFRA' : 'Obazoa',
                'COWC' : 'Obazoa',
                'MBRE' : 'Obazoa',
                'SROS' : 'Obazoa',
                'TADH' : 'Obazoa',
                'AQUE' : 'Obazoa',
                'ADIG' : 'Obazoa',
                'NVEC' : 'Obazoa',
                'HVUL' : 'Obazoa',
                'TKIT' : 'Obazoa',
                'MLEI' : 'Obazoa',
                'AVAG' : 'Obazoa',
                'EMUL' : 'Obazoa',
                'TSOL' : 'Obazoa',
                'TASI' : 'Obazoa',
                'GSAL' : 'Obazoa',
                'SMED' : 'Obazoa',
                'SMAN' : 'Obazoa',
                'PCAU' : 'Obazoa',
                'SRAT' : 'Obazoa',
                'SBAT' : 'Obazoa',
                'OVOL' : 'Obazoa',
                'PPAC' : 'Obazoa',
                'BMAL' : 'Obazoa',
                'CELE' : 'Obazoa',
                'TSPI' : 'Obazoa',
                'BANT' : 'Obazoa',
                'DMEL' : 'Obazoa',
                'AGAM' : 'Obazoa',
                'GMOR' : 'Obazoa',
                'BMOR' : 'Obazoa',
                'APIS' : 'Obazoa',
                'TCAS' : 'Obazoa',
                'RPRO' : 'Obazoa',
                'ZNEV' : 'Obazoa',
                'PHUM' : 'Obazoa',
                'ISCA' : 'Obazoa',
                'SMIM' : 'Obazoa',
                'SMAR' : 'Obazoa',
                'LSAL' : 'Obazoa',
                'DMAG' : 'Obazoa',
                'DPUL' : 'Obazoa',
                'LANA' : 'Obazoa',
                'HAZT' : 'Obazoa',
                'LPOL' : 'Obazoa',
                'RVAR' : 'Obazoa',
                'CGIG' : 'Obazoa',
                'LGIG' : 'Obazoa',
                'BGLA' : 'Obazoa',
                'OBIM' : 'Obazoa',
                'HROB' : 'Obazoa',
                'CTEL' : 'Obazoa',
                'ILIN' : 'Obazoa',
                'PFLA' : 'Obazoa',
                'SKOW' : 'Obazoa',
                'PMIN' : 'Obazoa',
                'SPUR' : 'Obazoa',
                'BBEL' : 'Obazoa',
                'BFLO' : 'Obazoa',
                'ODIO' : 'Obazoa',
                'CINT' : 'Obazoa',
                'BSCH' : 'Obazoa',
                'PETM' : 'Obazoa',
                'CMIL' : 'Obazoa',
                'DRER' : 'Obazoa',
                'TRUB' : 'Obazoa',
                'XTRO' : 'Obazoa',
                'APLA' : 'Obazoa',
                'PSIN' : 'Obazoa',
                'OANA' : 'Obazoa',
                'MMUS' : 'Obazoa',
                'HSAP' : 'Obazoa',
                'NUSP' : 'Obazoa',
                'FALB' : 'Obazoa',
                'RALL' : 'Obazoa',
                'MDAP' : 'Obazoa',
                'VCUL' : 'Obazoa',
                'EINT' : 'Obazoa',
                'PSPF' : 'Obazoa',
                'PSPE' : 'Obazoa',
                'NSPE' : 'Obazoa',
                'ASPE' : 'Obazoa',
                'OSPE' : 'Obazoa',
                'BDEN' : 'Obazoa',
                'SPUN' : 'Obazoa',
                'GPRO' : 'Obazoa',
                'AMAC' : 'Obazoa',
                'CANG' : 'Obazoa',
                'BBRI' : 'Obazoa',
                'BMER' : 'Obazoa',
                'CCOR' : 'Obazoa',
                'CREV' : 'Obazoa',
                'SMUC' : 'Obazoa',
                'LPEN' : 'Obazoa',
                'MPTE' : 'Obazoa',
                'RBRE' : 'Obazoa',
                'SPLU' : 'Obazoa',
                'MCIR' : 'Obazoa',
                'PBLA' : 'Obazoa',
                'LCOR' : 'Obazoa',
                'LTRA' : 'Obazoa',
                'MVER' : 'Obazoa',
                'MELO' : 'Obazoa',
                'RIRR' : 'Obazoa',
                'RSPE' : 'Obazoa',
                'UMAY' : 'Obazoa',
                'AING' : 'Obazoa',
                'MLAR' : 'Obazoa',
                'MOSM' : 'Obazoa',
                'ABIS' : 'Obazoa',
                'RSOL' : 'Obazoa',
                'CNEO' : 'Obazoa',
                'SCOM' : 'Obazoa',
                'SPOM' : 'Obazoa',
                'PMUR' : 'Obazoa',
                'SCER' : 'Obazoa',
                'KLAC' : 'Obazoa',
                'YLIP' : 'Obazoa',
                'NCRA' : 'Obazoa',
                'FOXY' : 'Obazoa',
                'TMEL' : 'Obazoa',
                'BHOM' : 'RASH',
                'AKER' : 'RASH',
                'SAGG' : 'RASH',
                'ALIM' : 'RASH',
                'PINF' : 'RASH',
                'HPAR' : 'RASH',
                'PHAL' : 'RASH',
                'PULT' : 'RASH',
                'ALAI' : 'RASH',
                'AAST' : 'RASH',
                'SPAR' : 'RASH',
                'NGAD' : 'RASH',
                'AANO' : 'RASH',
                'ESIL' : 'RASH',
                'COKA' : 'RASH',
                'PTRI' : 'RASH',
                'TPSE' : 'RASH',
                'FCYL' : 'RASH',
                'PMAR' : 'RASH',
                'SMIC' : 'RASH',
                'SMIN' : 'RASH',
                'PFAL' : 'RASH',
                'BBIG' : 'RASH',
                'TANN' : 'RASH',
                'GNIP' : 'RASH',
                'CMUR' : 'RASH',
                'CPAI' : 'RASH',
                'TGON' : 'RASH',
                'EACE' : 'RASH',
                'CVEL' : 'RASH',
                'VBRA' : 'RASH',
                'PTET' : 'RASH',
                'TTHE' : 'RASH',
                'IMUL' : 'RASH',
                'PPER' : 'RASH',
                'SLEM' : 'RASH',
                'OTRI' : 'RASH',
                'SCOE' : 'RASH',
                'BNAT' : 'RASH',
                'PBRA' : 'RASH',
                'RFIL' : 'RASH',
                'MONO' : 'Metamonada',
                'TVAG' : 'Metamonada',
                'GINT' : 'Metamonada',
                'SSAL' : 'Metamonada',
                'LMAJ' : 'Discoba',
                'TBRU' : 'Discoba',
                'EGRA' : 'Discoba',
                'ADEA' : 'Discoba',
                'SCUL' : 'Discoba',
                'PHYT' : 'Discoba',
                'PERK' : 'Discoba',
                'BSAL' : 'Discoba',
                'NGRU' : 'Discoba',
                'GTHE' : 'Archaeplastida',
                'EHUX' : 'RASH',
                'CTOB' : 'RASH',
                'CPAR' : 'Archaeplastida',
                'GSUL' : 'Archaeplastida',
                'CMER' : 'Archaeplastida',
                'PPUR' : 'Archaeplastida',
                'CCRI' : 'Archaeplastida',
                'CVAR' : 'Archaeplastida',
                'CSUB' : 'Archaeplastida',
                'CREI' : 'Archaeplastida',
                'MNEG' : 'Archaeplastida',
                'VCAR' : 'Archaeplastida',
                'MSPE' : 'Archaeplastida',
                'OLUC' : 'Archaeplastida',
                'KFLA' : 'Archaeplastida',
                'MPOL' : 'Archaeplastida',
                'PPAT' : 'Archaeplastida',
                'SFAL' : 'Archaeplastida',
                'SMOE' : 'Archaeplastida',
                'PABI' : 'Archaeplastida',
                'ATRI' : 'Archaeplastida',
                'OSAT' : 'Archaeplastida',
                'NNUC' : 'Archaeplastida',
                'ATHA' : 'Archaeplastida',
                'ACOE' : 'Archaeplastida'}

supergroups5 = {}
supergroups2 = {}
supergroups4 = {}
for sp, group in supergroups6.items():
    if group in ('Discoba', 'Metamonada'):
        supergroups5[sp] = 'Excavata'
        supergroups2[sp] = 'Diphoda'
        supergroups4[sp] = group
    elif group in ('Amoebozoa', 'Obazoa'):
        supergroups5[sp] = group
        supergroups2[sp] = 'Opimoda'
        supergroups4[sp] = 'Amorphea'
    else:
        supergroups5[sp] = group
        supergroups2[sp] = 'Diphoda'
        supergroups4[sp] = 'Diaphoretickes'

def get_root_daughters(root_name, supergroups4):
    """Classify the eukaryotic species according to the root. Recognised root names:
    - Opimoda-Diphoda (or OD)
    - Excavata (or EF)
    - Diaphoretickes (or Dia)
    - Metamonada (or M)
    - Discoba (or Dis)
    - ADis-DiaM
    - AM-DiaDis"""
    if root_name in ('Opimoda-Diphoda', 'OD'):
        root_daughters = {sp: ("Opimoda" if group == "Amorphea" else "Diphoda") for sp, group in supergroups4.items()}
    elif root_name in ('Excavata', 'EF'):
        root_daughters = {sp: ("Excavata" if group in ("Discoba", "Metamonada") else "non-Excavata") for sp, group in supergroups4.items()}
    elif root_name in ('Diaphoretickes', 'Dia'):
        root_daughters = {sp: (group if group == "Diaphoretickes" else "non-Diaphoretickes") for sp, group in supergroups4.items()}
    elif root_name in ('Metamonada', 'M'):
        root_daughters = {sp: (group if group == "Metamonada" else "non-Metamonada") for sp, group in supergroups4.items()}
    elif root_name in ('Discoba', 'Dis'):
        root_daughters = {sp: (group if group == "Discoba" else "non-Discoba") for sp, group in supergroups4.items()}
    elif root_name == "ADis-DiaM":
        root_daughters = {sp: ("ADis" if group in ("Amorphea", "Discoba") else "DiaM") for sp, group in supergroups4.items()}
    elif root_name == "AM-DiaDis":
        root_daughters = {sp: ("AM" if group in ("Amorphea", "Metamonada") else "DiaDis") for sp, group in supergroups4.items()}
    else:
        sys.stderr.write(f"Error: {root_name} not a valid root name.\n")
        root_daughters = None
    return root_daughters

def seq_info(species, eukarya_path = '/home/julian/julian2/snel-clan-genomes/eukarya_new'):
    seq_species_info = {}
    for spec in species:
        try:
            with open(f'{eukarya_path}/data_set/proteomes_metadata/{spec}.metadata.txt') as info_file:
                seq_species_info[spec] = {}
                info_file.readline()
                for line in info_file:
                    line = line.rstrip()
                    fields = line.split('\t')
                    if fields[7] != '1': # Not longest transcript
                        continue
                    seqid = fields[0]
                    gene_symbol = fields[2]
                    description = fields[4]
                    seq_species_info[spec][seqid] = (gene_symbol, description)
        except IOError:
            sys.stderr.write(f'Error: metadata file for {spec} not found.\n')
            return False
    return seq_species_info
