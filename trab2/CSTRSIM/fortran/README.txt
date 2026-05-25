# Necessary software: Fortran compiler, gnuplot


LINUX
# if not installed in Linux as administrator execute:
# sudo apt-get install gfortran
# sudo apt-get install gnuplot

WINDOWS
Free Fortran compiler: Force Fortran http://www.lepsch.com/downloads/Force3beta3Setup.exe
Gnuplot: https://sourceforge.net/projects/gnuplot/files/gnuplot/


# Compile the Fortran source
gfortran -O3 -o cstr cstr.f


# Execution produces data file 'X.csv'


# Execute script for fault 12
cstr < ./cfg/fault_12_abnormal_feed_flowrate.txt

# Visualize 14 variables and 4 constraints with 'gnuplot'
gnuplot gnuplot-data.gnu
# ... or with Python
python3 processX.py

#=================================

# Execute script for fault 17
cstr < ./cfg/fault_17_abnormal_jacket_effluent_pressure.txt

# Visualize 14 variables and 4 constraints with 'gnuplot'
gnuplot gnuplot-data.gnu
# ... or with Python
python3 processX.py output/X.csv
