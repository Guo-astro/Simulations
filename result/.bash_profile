# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/.local/bin:$HOME/bin


export PAT
alias src_bash="source /home/uchu/yansong.guo/.bash_profile"
alias ll="ls -al --block-size=M"
alias vl="vi stdout.log"
alias ebash="vi /home/uchu/yansong.guo/.bash_profile"
alias sbash="less /home/uchu/yansong.guo/.bash_profile"
alias cd_src="cd /xc/home/yansong.guo/FDPS_repo/FDPS-master/src"
alias modswp="module swap PrgEnv-cray PrgEnv-gnu"
alias cd_relax="cd /xc/home/yansong.guo/FDPS_repo/FDPS-master/sample/"
alias cd_init="cd /xc/home/yansong.guo/FDPS_repo/FDPS-master/POLYTROPE_TESTS/INIT_CONDS/"
alias cd_dens_relax="cd /xc/home/yansong.guo/FDPS-master/apps/HVCC_DENS_RELAX/"
alias cd_force_relax="cd /xc/home/yansong.guo/FDPS-master/apps/HVCC_FORCE_RELAX/"
alias cd_static_nonadiabatic="cd /xc/home/yansong.guo/FDPS-master/apps/HYDROSTATIC_NONADIABATIC_GSPH/"
alias cd_hvcc_chem="cd /xc/home/yansong.guo/FDPS-master/apps/HVCC_IMBH_NONADIABATIC_CHEM/"
alias cd_mhd_shocktube="cd /xc/home/yansong.guo/FDPS-master/apps/GSPMHD"
alias cp_result80="cp result/0080.dat /home/uchu/yansong.guo/src"
alias cd_result_relax="cd /xc/home/yansong.guo/FDPS_repo/FDPS-master/POLYTROPE_TESTS/BINOUT_FORCE_RELAX/"
alias cd_result_hvcc_bin="cd /xc/home/yansong.guo/FDPS_repo/FDPS-master/CO_0.2_0.4/bin"
alias cd_result_hvcc_vtk="cd /xc/home/yansong.guo/FDPS_repo/FDPS-master/CO_0.2_0.4/vtk"
alias cd_result_hvcc_ascii="cd /xc/home/yansong.guo/FDPS_repo/FDPS-master/CO_0.2_0.4/ascii"
alias cllog="rm -f stderr.log stdout.log"
alias mystat="squeue -u yansong.guo"
alias putjobL="sbatch runMPI_L.sh"
alias putjobM="sbatch runMPI_M.sh"
alias putjobS="sbatch runMPI_S.sh"
alias copyscp="scp -r  xcfront.yukawa.kyoto-u.ac.jp:/xc/home/yansong.guo/FDPS_repo/FDPS-master/CO_0.2_0.4/vtk/HVCC_H2_9_MBH5_imp_2.50/HVCC_chem_00[0-2][0-9].vtk /home/uchu/yansong.guo/FDPS_repo/CO_0.2_0.4/vtk"
alias revcopy="scp -r /home/uchu/yansong.guo/FDPS_repo/CO_0.2_0.4/vtk  xcfront.yukawa.kyoto-u.ac.jp:/xc/home/yansong.guo/FDPS_repo/FDPS-master/CO_0.2_0.4/vtk/HVCC_H2_9_MBH5_imp_2.50/HVCC_chem_00[0-2][0-9].vtk"
alias initcondcopy=" scp  /home/uchu/yansong.guo/FDPS_repo/0091.dat  xcfront.yukawa.kyoto-u.ac.jp:/xc/home/yansong.guo/FDPS_repo/FDPS-master/POLYTROPE_TESTS/INIT_CONDS/"
alias nochemcopy="scp -r  xcfront.yukawa.kyoto-u.ac.jp:/xc/home/yansong.guo/FDPS_repo/FDPS-master/CO_0.2_0.4/vtk/HVCC_MBH5_imp_2.50/HVCC_nochem_0000.vtk /home/uchu/yansong.guo/FDPS_repo/CO_0.2_0.4/vtk"
