set title "CSTR"
set xlabel "Time [min x sample interval]"
set ylabel "Variables"
input='output/cstr18.txt'
# firstline	It seems that gnuplot already ignores '#' lines
fl=1
set grid
set style data lines
set xrange [0:100]
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
plot	input every ::fl  using 0:($1/20 +0) title "c_{A0}",\
	input every ::fl  using 0:($2/0.25 +1)  title "Q_{1}",\
	input every ::fl  using 0:($3/30 +2)  title "T_{1}",\
	input every ::fl  using 0:($4/2 +3)  title "L",\
	input every ::fl  using 0:($5/2.8 +4)  title "c_{A}",\
	input every ::fl  using 0:($6/17.1 +5)  title "c_{B}",\
	input every ::fl  using 0:($7/80 +6)  title "T_{2}",\
	input every ::fl  using 0:($8/0.9 +7)  title "Q_{5}",\
	input every ::fl  using 0:($9/0.25 +8)  title "Q_{4}",\
	input every ::fl  using 0:($10/20 +9)  title "T_{3}",\
	input every ::fl  using 0:($11/10 +10)  title "h_{7}",\
	input every ::fl  using 0:($12/0.1 +11)  title "m_{1}",\
	input every ::fl  using 0:($13/0.6 +12)  title "m_{2}",\
	input every ::fl  using 0:($14/0.9 +13)  title "u_{2}",\
	input every ::fl  using 0:($15+1+ 14)  title "INV CON",\
	input every ::fl  using 0:(($16+0.15)+1+ 15)  title "MOLBAL CON",\
	input every ::fl  using 0:($17+1+ 16)  title "CW HL CON",\
	input every ::fl  using 0:($18+1+ 17)  title "EFFL HL CON"

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
