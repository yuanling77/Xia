import os

def run_simulation(replica_num):
    commands = [
        'pmemd -O -i min_steep.in -o min_steep.out -p Mpro-4wi_solvated.prmtop -c Mpro-4wi_solvated.inpcrd -r Mpro-4wi_min_steep.rst -ref Mpro-4wi_solvated.inpcrd -x Mpro-4wi_min_steep.mdcrd',
        'pmemd.cuda -O -i min_gradient.in -o min_gradient.out -p Mpro-4wi_solvated.prmtop -c Mpro-4wi_min_steep.rst -r Mpro-4wi_min_gradient.rst -x Mpro-4wi_min_gradient.mdcrd -ref Mpro-4wi_min_steep.rst',
        'pmemd.cuda -O -i nvt_5.0.in -o nvt_5.0.out -p Mpro-4wi_solvated.prmtop -c Mpro-4wi_min_steep.rst -r Mpro-4wi_nvt_5.0.rst -x Mpro-4wi_nvt_5.0.mdcrd -ref Mpro-4wi_min_steep.rst',
        'pmemd.cuda -O -i nvt_2.0.in -o nvt_2.0.out -p Mpro-4wi_solvated.prmtop -c Mpro-4wi_min_steep.rst -r Mpro-4wi_nvt_2.0.rst -x Mpro-4wi_nvt_2.0.mdcrd -ref Mpro-4wi_nvt_5.0.rst',
        'pmemd.cuda -O -i npt_2.0.in -o npt_2.0.out -p Mpro-4wi_solvated.prmtop -c Mpro-4wi_nvt_2.0.rst -r Mpro-4wi_npt_2.0.rst -x Mpro-4wi_npt_2.0.mdcrd -ref Mpro-4wi_nvt_2.0.rst',
        'pmemd.cuda -O -i npt_0.5.in -o npt_0.5.in.out -p Mpro-4wi_solvated.prmtop -c Mpro-4wi_npt_2.0.rst -r Mpro-4wi_npt_0.5.rst -x Mpro-4wi_npt_0.5.mdcrd -ref Mpro-4wi_npt_2.0.rst',
        'nohup pmemd.cuda -O -i md.in -o md.out -p Mpro-4wi_solvated.prmtop -c Mpro-4wi_npt_0.5.rst -r Mpro-4wi_md.rst -x Mpro-4wi_md.mdcrd -ref Mpro-4wi_npt_0.5.rst &'
    ]
    
    os.chdir(f'replica{replica_num}')
    for cmd in commands:
        os.system(cmd)
    os.chdir('..')

for i in range(1, 11):
    run_simulation(i)

 
 
