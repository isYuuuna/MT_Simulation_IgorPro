#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

function width(tott,maxstep,num,N)
variable tott,maxstep,num,N        // import variable for kpsim_rep function
wave xp,tp 
variable maxit = round(tott/maxstep) + 1  // max possible steps
variable t
variable i,j,k

make/o/n = (maxit) Mwidth // 1-D for storing width at each time point
make/o/n = (N) curr  // 1-D for storing current position of each MT

kpsim_rep(tott,maxstep,num,N) // run kpsim_rep first to get xp,tp

t = DimSize(xp,0)
for (i = 0; i < t; i ++) // for each time
	for (j = 0; j < N; j ++) // for each MT
		curr[j] = xp[i][j]
	endfor 
	Mwidth[i] = WaveMax(curr) - WaveMin(curr) // get the width
endfor 
	for(j = t; j < maxit; j ++)   // delete unfilled entries
		Mwidth[j] = nan
	endfor 

DoWindow /K MT_width
Display /N = MT_width

appendToGraph Mwidth vs tp[][0] // graph width vs acutal time
ModifyGraph rgb(Mwidth) = (128,0,0), lsize(Mwidth) = 2
for (i = 0; i < N; i ++) // compare width with xp vs tp 
	appendToGraph xp[][i] vs tp[][i]
	legend/W = MT_width "\s(xp#" + num2str(i) + ") MT" + num2str(i)
endfor

Legend

Label left "Width"
Label bottom "Time(s)"
end function



function posiForce(tott,maxstep,num,N,minF,maxF)
variable tott,maxstep,num,N,minF,maxF // minF,maxF: allow user input for force range imposed on the MT
wave xp,w
variable t // obtain the time length for the simulation
variable MT // obtain the number of MTs
variable rowSum // count the row sum of xp
variable i,j,k,z // counter 
variable maxit = round(tott/maxstep) + 1  
variable Mcol
String name = "force ="
variable together 
  
make/o/n = (maxit,(maxF - minF)) M // the column of the wave records mean position at each time step under a specific force
make/o/n = (N) xpc

for (i = 0; i < (maxF-minF); i ++) 
	w[10] = i + minF
	//print(i)
	//print(w[10]) // uncomment the code to check the value of i and w[10]
	kpsim_rep(tott,maxstep,num,N) // run the simulation under different force to get different xp 
	t = DimSize(xp,0)
	MT = DimSize(xp,1)    
	for(j = 0; j < t; j ++)
		rowSum = 0
		for (k = 0; k < MT; k ++)
			rowSum = rowSum + xp[j][k]
			xpc[k] = xp[j][k]
		endfor
    //WaveMax(xpc) = max position 
	M[j][i] = rowSum/N
	endfor
	for(j = t; j < maxit; j++)
		M[j][i] = nan
	endfor
	//deletepoints t,(maxit - t),M[][i]
	//
endfor
	
	Mcol = DimSize(M,1)
	
	DoWindow /K Yuna_test
	Display /N =Yuna_test
	
	for(i = 0; i < Mcol; i ++) // plot the graph of mean position vs time under different force
	    //together = name + num2str(i+minF)
		appendToGraph/C =(500*i, 15500*i,0), M[][i]
		legend
		//legend \s(copper) copper
	endfor
	
	ModifyGraph mirror=2
    ModifyGraph fStyle=1
    ModifyGraph axThick=1.5
    Label left "mean position"
    Label bottom "time"

end function

function posiStiff(tott,N,maxstep,minF,maxF)
variable tott,N,maxstep,minF,maxF // minF,maxF: allow user input for force range imposed on the MT
wave xp,w
variable t = DimSize(xp,0) // obtain the time length for the simulation
variable MT = DimSize(xp,1) // obtain the number of MTs
variable rowSum // count the row sum of xp
variable i,j,k // counter 
  
make/o/n = (t,(minF - maxF + 1)) M // the column of the wave records mean position at each time step under a specific force

for (i = minF; i < maxF; i ++)
	w[10] = i
	kpsim(tott,maxstep,N)     // run the simulation under different force to get different xp 
	for(j = 0; i < t; i ++)
		rowSum = 0
		for (k = 0; j < MT; j ++)
			rowSum = rowSum + xp[j][k]
		endfor
	M[j][i] = rowSum/N
	endfor
endfor

end function
