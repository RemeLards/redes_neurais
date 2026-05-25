set title "CSTR"
set xlabel "Time [min x sample interval]"
set ylabel "Variables"
input='./output/X.csv'
# firstline
fl=2
set grid
set style data lines
set yrange [-1:20]
set key right top outside
#set autoscale y


# http://xmodulo.com/how-to-plot-using-specific-rows-of-data-file-with-gnuplot.html
#
# plot "my.dat" every A:B:C:D:E:F
# 
# A: line increment
# B: data block increment
# C: The first line
# D: The first data block
# E: The last line
# F: The last data block
# To plot the data from line 10 to line 100: plot "my.dat" every ::10::100
#
# Last four values are constraints with nominal value zero
plot	input every ::fl  using 0:($1/20+0) title "c_{A0}",\
	input every ::fl  using 0:($2*4+1)  title "Q_{1}",\
	input every ::fl  using 0:($3/30+2)  title "T_{1}",\
	input every ::fl  using 0:($4/2+3)  title "L",\
	input every ::fl  using 0:($5/3+4)  title "c_{A}",\
	input every ::fl  using 0:($6/17+5)  title "c_{B}",\
	input every ::fl  using 0:($7/80+6)  title "T_{2}",\
	input every ::fl  using 0:($8*1+7)  title "Q_{5}",\
	input every ::fl  using 0:($9*4+8)  title "Q_{4}",\
	input every ::fl  using 0:($10/20+9)  title "T_{3}",\
	input every ::fl  using 0:($11/60000 + 10)  title "PCW",\
	input every ::fl  using 0:($12/25+11)  title "CNT_{1}",\
	input every ::fl  using 0:($13/40+12)  title "CNT_{3}=SP_Q5",\
	input every ::fl  using 0:($14*1+13)  title "CNT_{2}",\
	input every ::fl  using 0:($15*1+15)  title "INV CON",\
	input every ::fl  using 0:($16*1+16)  title "CW DP CON",\
	input every ::fl  using 0:($17*1+17)  title "EFFL DP CON",\
	input every ::fl  using 0:($18*1+18)  title "MOLBAL CON"

print "Hit <Enter> to continue..."
pause -1

set terminal png 
print "Dumping PNG image file..."
set output "./output/Xgnuplot.png"
replot


# 01 FEED_CONCENTRATION_SENSOR       1989 6 1 0 0 0       19.6676KMOL/M3
# 02 FEED_FLOWRATE_SENSOR            1989 6 1 0 0 0        0.2502M3/MIN 
# 03 FEED_TEMPERATURE_SENSOR         1989 6 1 0 0 0       29.8533C      
# 04 REACTOR_LEVEL_SENSOR            1989 6 1 0 0 0        2.0057M      
# 05 CONCENTRATION_A_SENSOR          1989 6 1 0 0 0        2.8862KMOL/M3
# 06 CONCENTRATION_B_SENSOR          1989 6 1 0 0 0       17.1528KMOL/M3
# 07 REACTOR_TEMPERATURE_SENSOR      1989 6 1 0 0 0       79.9374C      
# 08 COOLING_WATER_FLOWRATE_SENSOR   1989 6 1 0 0 0        0.9213M3/MIN 
# 09 PRODUCT_FLOWRATE_SENSOR         1989 6 1 0 0 0        0.2526M3/MIN 
# 10 COOLING_WATER_TEMPERATURE_SENSOR1989 6 1 0 0 0       19.9116C      
# 11 COOLING_WATER_PRESSURE_SENSOR   1989 6 1 0 0 0    55893.9436KG/M2  
# 12 LEVEL_CONTROLLER_OUTPUT_SIGNAL  1989 6 1 0 0 0       25.3117% OPEN 
# 13 CW_FLOW_CONTROLLER_OUTPUT_SIGNAL1989 6 1 0 0 0       40.6786% OPEN 
# 14 COOLING_WATER_SETPOINT          1989 6 1 0 0 0        0.9000M3/MIN 
# 15 INVENTORY CONSTRAINT            1989 6 1 0 0 0        0.0086M3     
# 16 CW_PRESSURE_DROP_CONSTRAINT     1989 6 1 0 0 0        0.0250M3/MIN 
# 17 EFFL_PRESSURE_DROP_CONSTRAINT   1989 6 1 0 0 0        0.0021M3/MIN 
# 18 MOL_BALANCE_CONSTRAINT          1989 6 1 0 0 0        0.3600KMOL   
