set title "CSTR"
set xlabel "Time [min x sample interval]"
set ylabel "Variables"
input='output/cstr18.txt'
# firstline	It seems that gnuplot already ignores '#' lines
fl=1
set grid
set style data lines
set xrange [0:100]
set yrange [-1:5]
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
plot	input every ::fl  using 0:($15+1+ 0)  title "INV CON",\
	input every ::fl  using 0:(($16+0.15)+1+ 1)  title "MOLBAL CON",\
	input every ::fl  using 0:($17+1+ 2)  title "CW HL CON",\
	input every ::fl  using 0:($18+1+ 3)  title "EFFL HL CON"

print "Hit <Enter> to continue..."
pause -1


#    1. Feed concentration
#    2. Feed flowrate
#    3. Feed temperature
#    4. Reactor level
#    5. Product A concentration
#    6. Product B concentration
#    7. Reactor temperature
#    8. Coolant flowrate
#    9. Product flowrate
#    10. Coolant inlet temperature
#    11. Coolant inlet head
#    12. Level controller output
#    13. Coolant controller output
#    14. Coolant setpoint
#    15. Inventory
#    16. Mol balance
#    17. Cooling water head loss
#    18. Effluent head loss
#    19. Class
