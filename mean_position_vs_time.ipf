#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

function decimate(input, factor)
wave input
variable factor
make/o/n = (numpnts(input)/factor) decimated = nan

decimated = input[p*factor]
end function

// This tdetector function is based on the MATLAB code in 'Molecular counting by photobleaching...'
// Because this algorithm declares way too many steps, we could change the z-score multiplier as a way
// to decrease the number of declared steps. 
function tdetector(X, VO)
wave X // vector of a piecewise constant function hidden in white noise
variable VO
variable Lo,i,platsCount,step_index,status,foundCount, j
variable maxit = 128

Lo = numpnts(X) - 1
varSect(X, maxit)
wave vx
make/o/n = (maxit) Ywave // the final output wave
make/o/n = (Lo + 1) multTabWave 

wave multTab1,multTabx,multTaby
for(i = 0; i < Lo + 1; i ++)
	multTabWave = interp(i,multTabx,multTaby)
endfor

make/o/n = (maxit,2) plats_array = nan
make/o/n = (maxit) found = nan
make/o/n = (2) Bound
make/o/n = (maxit) tool
make/o/n = (maxit) tool1

plats_array[0][0] = 0
plats_array[0][1] = Lo
platsCount++ // platsCount = 1

do
 Bound[0] = plats_array[platsCount - 1][0]
 Bound[1] = plats_array[platsCount - 1][1]
 
 Duplicate/O/R=[Bound[0], Bound[1]] X, tool
 Duplicate/O/R=[Bound[0], Bound[1]] vx, tool1
 [step_index,status] = detectStep2(tool, Bound[0], tool1, multTabWave)
 
 if (status == 1)
 	found[foundCount] = step_index
 	foundCount++
 	plats_array[platsCount][0] = plats_array[platsCount - 1][0]
 	plats_array[platsCount][1] = step_index - 1
 	platsCount++
 	plats_array[platsCount][0] = step_index
 	plats_array[platsCount][1] = plats_array[platsCount - 2][1]
 	platsCount++
 	plats_array[platsCount - 3, ] = nan
 	platsCount--	
 elseif (status == -1)
 	plats_array[platsCount - 1, ] = nan 
 	platsCount--	
 endif
while(platsCount != 0)

make/o/n = (1) one
one = 1
make/o/n = (1) Lo_1
Lo_1 = Lo
Concatenate/O {found, one, Lo_1}, found
Sort found,found

checkSteps(found, X, vx, multTabWave, VO)
wave checked
Lo_1 = Lo_1 + 1
Concatenate/O {one, checked,Lo_1}, checked

for(i = 0; i < numpnts(checked) - 1; i ++)
	for(j = checked[i]; j < checked[i + 1] - 1; j ++)
		Ywave[j] = mean(X, checked[i], checked[i + 1])
	endfor
endfor

make/o/n = (numpnts(checked)-3) step_sizes = 0
for(i = 1; i < numpnts(checked); i ++)
	step_sizes[i - 1] = Ywave[checked[i]] - Ywave[checked[i] - 1]
endfor

end function

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function getSig(Xs, expnt, maxit) // output: SIG
Wave Xs
variable expnt,maxit
variable sigmaC, sigmaN, true, SIG = 0
variable i,j,k
make/o/n = (maxit) diff1
make/o/n = (maxit) diff2

for (i = 0; i < maxit - 1; i ++)
	diff1 = Xs[i+1] - Xs[i]
endfor 
diff2 = diff1 * diff1

do
sigmaC = sqrt(mean(diff2)/2)
for (i = 0; i < maxit - 1; i ++)

	if (diff2[i] > 3*(2^0.5)*sigmaC)
		diff2[i] = nan
	endif
	if (diff1[i] > 3*(2^0.5)*sigmaC)
		diff1[i] = nan
	endif
	
	sigmaN = sqrt(mean(diff2)/2)
	
	if (sigmaN == sigmaC)
		break
	endif
endfor 
while (true == 0)

SIG = sigmaN * 1.015
SIG = SIG^expnt

return SIG
end function

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function getSigLoop(Xs, expnt, sd2, ld2, maxit) // output: SIG
Wave Xs
variable expnt, sd2, ld2, maxit
variable sigmaC, sigmaN, true, sts, lts, STOP, SIG = 0
variable i,j,k
make/o/n = (maxit) diff1
make/o/n = (maxit) diff2
make/o/n = (maxit) icurrpeaks
make/o/n = (maxit) currpeaks

getDiff(Xs)
wave diff
diff1 = diff 
diff2 = diff1 * diff1

do
sigmaC = sqrt(mean(diff2)/2)
	SIG = ((sd2-sts)/(ld2-lts)/2)^0.5	
	sd2 = sd2-sts
    ld2 = ld2-lts
    
    for (i = 0; i < numpnts(diff1); i ++)
    	if (abs(diff1[i]) > (3*(2^.5))*SIG)
    		icurrpeaks[i] = i
    	endif
    endfor
    
    for (i = 0; i < numpnts(icurrpeaks); i ++)
    	currpeaks[i] = diff2[icurrpeaks[i]]
    endfor
    
    for(i = 0; i < numpnts(icurrpeaks); i ++)
    	diff1[icurrpeaks[i]] = 0
    endfor 
    
    sts = sum(currpeaks)
    lts = numpnts(currpeaks)
    
    if (lts == 0)
    	break
    endif	
while (STOP == 0)

SIG = sigmaN * 1.015
SIG = SIG^expnt

return SIG
end function
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function varSect(X, maxit) // output: vx
wave X
variable maxit
variable i,j,Lo,step_index,status,foundCount,platsCount

foundCount = 0
make/o/n = (maxit,2) plats_array = nan
make/o/n = (maxit) found = nan
make/o/n = (2) Bound
make/o/n = (maxit) tool
make/o/n = (maxit) vx

Lo = numpnts(X) - 1
plats_array[0][0] = 0
plats_array[0][1] = Lo
platsCount++ // platsCount = 1

do
 Bound[0] = plats_array[platsCount - 1][0]
 Bound[1] = plats_array[platsCount - 1][1]
 
 Duplicate/O/R=[Bound[0], Bound[1]] X, tool
 [step_index,status] = detectVars(tool,Bound[0],maxit)
 
 if (status == 1)
 	found[foundCount] = step_index
 	foundCount++
 	plats_array[platsCount][0] = plats_array[platsCount - 1][0]
 	plats_array[platsCount][1] = step_index - 1
 	platsCount++
 	plats_array[platsCount][0] = step_index
 	plats_array[platsCount][1] = plats_array[platsCount - 2][1]
 	platsCount++
 	plats_array[platsCount - 3, ] = nan
 	platsCount--	
 elseif (status == -1)
 	plats_array[platsCount - 1, ] = nan 
 	platsCount--	
 endif
while(platsCount != 0)

make/o/n = (1) one
one = 1
make/o/n = (1) Lo_1
Lo_1 = Lo
Concatenate/O {found, one, Lo_1}, found
Sort found,found

checkVars(found, X, maxit)
wave checked
Lo_1 = Lo + 1
Concatenate/O {one, checked, Lo_1}, checked

for (i = 0; i < numpnts(checked); i ++)
	for(j = checked[i]; j < checked[i + 1]; j ++)
		Duplicate/O/R=[checked[i], checked[i + 1] - 1] X, help
		vx[j] = getSig(help, 2, maxit)
	endfor
endfor 

end function
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [Variable mxiNum, Variable status] detectVars (wave Xs, variable i_1, variable maxit)
variable i,j,k
variable L , Asd2, Bsd2, SIG, VA, VB, DOV, LA, LB, sigma_squared
make/o/n = (maxit) mxi

L = numpnts(Xs) - 1
status = -1 
mxi = 0

if (L >= 22)
	getDiff(Xs)
	wave diff
	make/o/n = (numpnts(diff)) d2
	d2 = diff^2
	
	SIG = getSig(Xs,1,maxit)

	Asd2 = sum(d2,0,9)
	Bsd2 = sum(d2,10,numpnts(d2))
	make/o/n = (L-2) RVD = 0
	
	for (i = 10; i < L - 9; i ++)
		Asd2 = Asd2 + d2[i-1];
        Bsd2 = Bsd2 - d2[i];
        
        Duplicate/O/R=[0,i] Xs,A
		Duplicate/O/R=[i+1,inf] Xs,B
		
		VA = getSigLoop(A,2,Asd2,i-1,maxit)
        VB = getSigLoop(B,2,Bsd2,L-i-1,maxit)       
        DOV = VA - VB;
        
        LA = i
        LB = L - i
        sigma_squared = ( ((LA^2 + LA -3)/((LA-1)^2)) + ((LB^2 + LB -3)/((LB-1)^2)) - 2 )*SIG^4;
        RVD[i] = DOV/(sigma_squared^.5)/3;	
	endfor
	
	make/o/n = (L-2) absRVD
	absRVD = abs(RVD)
	Extract/INDX absRVD, mxi, absRVD == WaveMax(absRVD)
    mxiNum = WaveMax(mxi) 

	 if (abs(RVD(mxiNum)) > 1)
        status = 1;
        mxiNum = mxiNum + i_1;   
    endif
	
endif 

return [mxiNum, status]
end function

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function getDiff(X)
wave X
variable i, L

L = numpnts(X)
make/o/n = (L) diff

for(i = 0; i < L - 1; i ++)
	diff[i] = X[i + 1] - X[i]
endfor
end function

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function checkVars(wave found, wave rx, variable maxit)
variable ii, cc, endW,step_index, status

found[numpnts(found) - 1] = found[numpnts(found) - 1] + 1
make/o/n = (120001) checked = nan
make/o/n = (120001) rx
cc = 0
endW = 0

if(numpnts(found) - 1 == 2)
	endW = 1
endif

ii = 0
do
	ii++
	Duplicate/O/R=[found[ii - 1],found[ii + 1] - 1] found,rx
	[step_index, status] = detectVars(rx, found[ii - 1], maxit)

	if (status == 1)
		cc = cc + 1;
		checked[cc] = step_index
    else
        found[ii] = nan;
        ii = ii - 1;        
    endif
    
    if(ii + 1 == numpnts(found) - 1)
    	endW = 1
    endif
while(endW == 0)

end function

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function [variable mxiNum, variable status] detectStep2(wave Xs, variable i_1, wave vx, wave multTab)
variable L,SIG,mult,m1,m2,ii,DOM,sigma,sigma_vx, RMD_vx

L = numpnts(Xs) - 1 
status = -1
mxiNum = 0

if (L >= 2) 
    SIG = mean(vx)^.5;
    mult = multTab(L);
    
    make/o/n = (L) RMD 
    RMD = 0
    
    m1 = 0;
    m2 = sum(Xs);
    for (ii = 0; ii < L - 1; ii ++)
        m1 = m1 + Xs[ii]
        m2 = m2 - Xs[ii]
        DOM = m2/(L-ii) - m1/(ii)
        sigma = SIG*(1/ii + 1/(L-ii))^.5
        RMD[ii+1] = DOM/(sigma*mult);
    endfor
    
    make/o/n = (L) absRMD 
    absRMD = abs(RMD)
	Extract/INDX absRMD, mxi, absRMD == WaveMax(absRMD)
    mxiNum = WaveMax(mxi) 
endif 

sigma_vx = sum(vx,0,mxiNum)/(mxiNum - 1)^2 + sqrt((sum(vx,mxiNum,inf))/(L - mxiNum + 1)^2)
DOM = mean(Xs,0,mxiNum) - mean(Xs,mxiNum,inf)
RMD_vx = DOM/(mult*sigma_vx)

 status = -1;
 if (abs(RMD[mxiNum]) > 1 && abs(RMD_vx) > 1)
    status = 1;
 	mxi = mxi + i_1 - 1;   
 endif

end function
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function checkSteps (found, rx, noise_input, multTabWave, VO)
wave found, rx, multTabWave, noise_input
variable VO
variable cc, endW,step_index, status,ii

make/o/n = (1200001) vx = nan
vx = noise_input
found[numpnts(found) - 1] = found[numpnts(found) - 1] + 1
make/o/n = (1200001) checked = nan
cc = 0
endW = 0

if(numpnts(found) - 1 == 2)
	endW = 1
endif

ii = 0
do
	ii++
	Duplicate/O/R=[found[ii - 1],found[ii + 1] - 1] rx, tool1 
	Duplicate/O/R=[found[ii - 1],found[ii + 1] - 1] vx, tool2
	[step_index, status] = detectStep2(tool1, found[ii - 1], tool2, multTabWave)

	if (status == 1)
		cc = cc + 1;
		checked[cc] = step_index
    else
        found[ii] = nan;
        ii = ii - 1;        
    endif
    
    if(ii + 1 == numpnts(found) - 1)
    	endW = 1
    endif
while(endW == 0)

end function

//---------------------------------------------------------------------------------------------------------------------------------------------------------
function help1()
wave Std 
Std = 10000*Std
end function

function help()
variable noise,i
wave filter3
variable maxit = 1200001
make/o/n = (maxit) DOM // difference of variance
make/o/n = (maxit) subset1 = nan
make/o/n = (maxit) subset2 = nan
for (i = 0; i < maxit ; i++)
	Duplicate/O/R=[0,i-1] filter3, subset1
	Duplicate/O/R=[i,maxit - 1] filter3, subset2
	DOM[i] = mean(subset1) - mean(subset2)
endfor
end function

function cal()
wave DOV5000,Std
variable maxit = 1200001
make/o/n = (1200001) absDOV5000
variable i,noise
noise = 2.6180*10^(-12)
for (i = 0; i < maxit ; i++)
	Std[i] = (noise^2)*((i^2 + i - 3)/(i - 1)^2 + ((maxit - i)^2 + (maxit - i) - 3)/((maxit - i)^2 - 1)^2 - 2)	// recheck 
endfor
end function

function detector()
//    Tdetector1
// 1. Calculate the underlying white noise sigma^2.
// 2. Iterate through every possible way of splitting X into two sections, and calculate the DOM.
// 3. Calculate the significance of DOM by taking the absolute value of (multiplier * sd of its respective DOM distribution)
// 4. Only the most significance DOM results in a declared step.
// 5. Iterate the process until no new plateaus are declared. 

//    Tdetector2
// 1. Iterate through every possible way of splitting X into two sections, and calculate the DOV.
// 2. Calculate the distribution of DOV using the 'current' subset of the input vector.
// 3. Calculate the significance of DOV. ???
// 4. Only the most significance DOV results in a declared step.
// 5. Iterate the process until no new plateaus are declared. 
// 6. Repeat the DOM detection in Tdetector1 with adjusted distribution.
//wave KTv
wave filter3,tp1
variable maxit = 1200001
variable i,j,k,L,multipler,count,anySectionLeft
variable T,P
variable significance // only takes 0 or 1
make/o/n = (maxit) DOV // difference of variance
make/o/n = (maxit) absDOV
make/o/n = (maxit) DOM // difference of variance
make/o/n = (maxit) Std // standard deviation of corresponding DOV 
make/o/n = (maxit) subset1 = nan
make/o/n = (maxit) subset2 = nan
make/o/n = (maxit,2) step // time + position

L = maxit
T = 0
P = 0
count = 0
anySectionLeft = 1

significance = 0
variable noise // constant: this is the variance of underlyinng noise throughout the input vector
for (i = 0; i < L - 1; i++)
	noise += (filter3[i+1] - filter3[i])^2
endfor 
noise = noise/(2*(L - 1))


for (i = 0; i < L ; i++)
	Duplicate/O/R=[0,i-1] filter3, subset1
	Duplicate/O/R=[i,L - 1] filter3, subset2
	DOV[i] = variance(subset1) - variance(subset2)
	Std[i] = (noise^2)*((i^2 + i - 3)/(i - 1)^2 + ((L - i)^2 + (L - i) - 3)/((L - i)^2 - 1)^2 - 2)	
endfor
multipler = -sqrt(2)*InverseERFC(-0.95^(1/(L - 1)))
absDov = abs(DOV)
for (i = 0; i < L; i ++)
	if (absDOV[i] > multipler * Std[i])
		significance = 1
		if (absDOV[i] > P)
			P = absDOV[i]
			T = tp1[i]
		endif 
	endif
endfor
step[count][0] = T
step[count][1] = P
count++
anySectionLeft = 2*anySectionLeft
end function

//-----------------------------------------------------------------------------------------------------------------------------------------------------
function GaussianFilter(tott,maxstep)
// G(x) = 1/(sqrt{2pi}*sigma)*exp(-x^2/(2*sigma^2))
variable tott,maxstep
variable maxit = round(tott/maxstep) + 1
variable kernelWidth, radius, i, j, k, count, num, sigma,kernelSum,velocity 
wave KTv5

sigma = 500
radius = 6*sigma
kernelWidth = 2*radius + 1 // 50 seconds time window
num = 0
make/o/n = (kernelWidth + 1) kernel
make/o/n = (maxit) filter6

// create kernel
count = 0
for (j = - radius; j < radius + 1; j ++) // increment of 1 second
	kernel[count] = (1/(sqrt(2*pi)*sigma))*exp(-(j^2)/(2*sigma^2))
	kernelSum += kernel[count]
	count ++
endfor
kernel = kernel/kernelSum

for(i = 0; i < maxit; i ++) 
	count = 0
	for (j = i - radius; j < i + radius + 1; j ++) 
		if (j < 0 || j >= maxit)
			velocity = 0
		else 
			velocity = KTv5[j]
		endif 
		filter6[num] += kernel[count]* velocity
		count++
	endfor
	num++
endfor

end function
//-------------------------------------------------------------------------------------------------------------------------------------------------
function change()
wave pointPos,pointSize
variable i,j,k,q

for(i = 0; i < 42; i ++)
	 if(pointSize[i]/20 < 3)
	 	ModifyGraph msize(pointPos[i])= 3
	 else
	 	ModifyGraph msize(pointPos[i])= pointSize[i]/20
	 endif
endfor

end function
//---------------------------------------------------------------------------------------------------------------------------------------------------
function Fig2_KT(tott,maxstep,num,N)
variable tott, maxstep, num, N
variable maxit = round(tott/maxstep) + 1
wave w,o,xp,tp,vp,M,final
variable i,j,k,q
variable endT,startT,endI,startI,flag,counter,total,ap,sameBranch// if KTVswitch= 1, use kinechore veloctiy; otherwise, use mean MT bundle velocity

Make/o/n = (maxit) tp1,xp1 // 1-D wave tp
Make/o/n=(2) average = 0 // the first entry records average velocity; the second entry record waiting time for next switch
Make/o/n=(2) single = 0
Make/o/n=(1000,2) pointPos = nan
Make/o/n=(1000) pointSize = nan 

// set parameters
w[8] = 20// set the link stiffness to 20/5/7.5 pN/um
o[12] = 1 // turn on the viscous option
o[4] = 0 // set vary assembly velocity to 0
o[6] = 0 // set vary assembly velocity to 0
total = 0 // counter for pointPos and pointSize
 
for(ap = 0; ap < 1; ap ++) // make MTs to start at 2 different branches 
	print(ap)
	o[13] = ap // ap = 0: MTs are all disassembly; ap = 1: MTs are all assembly
	for (i = 5; i < N+1; i ++)
		print(i)
		w[10] = i // set the external force
		average = 0 // initialize 
		
		for(j = 0; j < num; j ++) // for each simulation
			single = 0
			counter = 0 // count the number of data points for each velocity; initialize
			sameBranch = 0 // 0 = not in the same branch		
			kpsim_initialize(tott,maxstep,N)
			kpsim(tott,maxstep,N)
			for(k = 0; k < maxit; k ++)
				xp1[k] = xp[k][0]
				tp1[k] = tp[k][0]
			endfor 
			getKTvelocity(tott,maxstep,N,xp1) // calculate the KTv,KTp based on xp and tp
			wave KTp,KTv
			getPeakV3(tott,maxstep,num,N,KTp)
			WaveStats final 
			
			flag = final[0][2] // 1 = assembly; -1 = disassembly
			
			if ((ap == 0 && flag == -1) || (ap == 1 && flag == 1)) 
				sameBranch = 1 // confirm that they are on the same branch
			endif 
			
			if(V_npnts/3 == 1) // V_npnts/3 
				if (sameBranch == 1) // if they are at the same branch 
					startT = final[0][0]
					findValue/T=0.002 /V=(startT) tp1
					startI = V_value // change the start time to corresponding index in tp1
					endI = maxit // used to be (maxit/2), which is an error. On line 75, we will cut the index in half. 
				else
					startI = 0
					findValue/T=0.002 /V=(final[0][0]) tp1
					endI = V_Value
				endif 
			else
				if (sameBranch == 1)
					startT = final[0][0] // start time
					endT = final[1][0] // end time
					findValue/T=0.002 /V=(startT) tp1
					startI = V_value // change the start time to corresponding index in tp1
					findValue/T=0.002 /V=(endT) tp1
					endI = V_value // change the end time to corresponding index in tp1 
				else
					startI = 0
					findValue/T=0.002 /V=(final[0][0]) tp1
			    	endI = V_Value
				endif 
			endif
			
			for(q = startI; q < round(endI/2); q ++) // only include the first half of the velocity 
				single[0] += KTv[q]
				single[1] += round(endT/2)-startT // total time it takes to switch
				counter++
			endfor
			print(single[1])
			print(single[0])
			if(counter == 0)
				single[0] = 0
			else
				single[0] = single[0]/counter
			endif
			print(single[0])
				
			average[0] += single[0]
			print(average[0])
			average[1] += (single[1]/(maxit/2))
		endfor
		
		average[0] = average[0]/num
		print(average[0])
		average[1] = average[1]/num
		print(average[1])
		pointPos[total][0] = average[0]// store the velocity
		pointPos[total][1] = i/N // store (External Force/Number of MTs)
		pointSize[total] = average[1]
		total++
	endfor 
endfor 

DoWindow /K Fig2_KT
Display /N = Fig2_KT
appendToGraph pointPos[][0] vs pointPos[][1] // make the basic graph 
ModifyGraph mode=3,marker=18
Label bottom "External Force/MT [pN]"
Label left "KT velocity[um/s]"
end function 

//-----------------------------------------------------------------------------------------------------------------------------------------------
// return the KTp velcity along each time point
function getKTvelocity(tott,maxstep,N,target)
variable N,tott,maxstep
wave target
variable maxit = round(tott/maxstep) + 1  
wave xp,tp,vp,force,w,tp1
variable i,j

make/o/n = (maxit) KTp  // wave that stores KT position at each time point (should be the same for all N MT-tips)
make/o/n = (maxit) KTv  // wave that stores KT velocity at each time point (should be the same for all N MT-tips)
make/o/n = (maxit) tp1

for(j = 0; j < maxit; j ++)
		tp1[j] = tp[j][0]
endfor 

for (i = 0; i < maxit; i ++)
	KTp[i] = target[i] + (force[i][0]/w[8])
endfor 

Differentiate KTp /X = tp1 /D = KTv
end function 
//--------------------------------------------------------------------------------------------------------------------------
function Fig2(tott,maxstep,num,N)
variable tott, maxstep, num, N
variable maxit = round(tott/maxstep) + 1
wave w,o,xp,tp,M
variable i,j,z,p,q,pos,neg,counter
string name
variable KTVswitch // if KTVswitch= 1, use kinechore veloctiy; if KTV = 0, use mean MT bundle velocity

w[8] = 20 // set the link stiffness to 20 pN/um
o[12] = 1 // turn on the viscous option
counter = 0
KTVswitch = 1 // indicate whether or not we want to use kinechore velocity in the following code
Make/o/n=(maxit) stable = nan
Make/o/n=(maxit) tpCopy
Make/o/n=(2) average = 0
Make/o/n=(1000,2) pointPos = nan
Make/o/n=(1000) pointSize = nan 

for(i = -N; i < N+1; i += 1)
	average = 0 // initialize
	pos = 0 // counter for the number of data points for positive velocity; initialize
	neg = 0 // counter for the number of data points for negative velocity; initialize
	w[10] = i // set the external force
	meansOfMeansMatrix(tott, maxstep,num,N)
	
	for(z = 0; z < maxit; z++) // make 1d wave for tp 
		tpCopy[z] = tp[z][0]
	endfor
	
	Differentiate M /X = tpCopy// get the velocity while set the xwave to tpCopy
	
	if (KTVswitch == 1)
		getKTvelocity(tott,maxstep,N,M)
		wave KTv,KTp
	endif 

	for(j = 20000; j < maxit; j++)
		if (KTVswitch == 1)
			stable[j - 20000] = KTv[j]
		else
			stable[j - 20000] = M[j]
		endif 
		if(stable[j - 20000] > 0) // count the positve and negative velocity
			average[0] += stable[j - 20000]
			pos++
		else
			average[1] += stable[j - 20000]
			neg++
		endif
	endfor
	
	if (pos == 0)
		average[0] = 0
	else
		average[0] = average[0]/pos // mean value for negative velocity 
	endif   
    if (neg == 0)
    	average[1] = 0
    else 
    	average[1] = average[1]/neg // mean value for positive velocity 
    endif 
	p = (neg/380000)*100 // the percentage of negative velocity 
	q = (pos/380000)*100 // the percentage of positve velocity 
	pointPos[counter][0] = average[0] // store the positive velocity
	pointPos[counter][1] = i/N // store (External Force/Number of MTs)
	pointSize[counter] = q
	pointPos[counter + 1][0] = average[1] // repeat the above process for negative velocity
	pointPos[counter + 1][1] = i/N 
	pointSize[counter + 1] = p
	counter = counter + 2 
endfor

DoWindow /K Fig2 
Display /N = Fig2
appendToGraph pointPos[][0] vs pointPos[][1] // make the basic graph 
ModifyGraph mode=3,marker=18

for(i = 0; i < counter; i ++)
	if(pointPos[i][0] != nan)
		if(pointSize[i][0] < 15)
			ModifyGraph msize(pointPos[i])= 1.5
		else
			ModifyGraph msize(pointPos[i])= pointSize[i]/10
		endif	
	endif	 
endfor

Label bottom "External Force/MT [pN]"
Label left "MT velocity[um/s]"

end function 
//-------------------------------------------------------------------------------------------------------------------------
function getHist(tott,maxstep,num,N)
variable tott, maxstep, num, N
variable maxit = round(tott/maxstep) + 1
wave w,o,xp,tp,M
variable i,j,z,p,q,pos,neg
string name

w[8] = 20 // set the link stiffness to 20 pN/um
o[12] = 1 // turn on the viscous option
Make/o/n=(maxit) stableM = nan
Make/o/n=(maxit) tpCopy
Make/o/n=(2) average
Make/o/n=(2) bin
bin[0] = -50
bin[1] = 50

for(i = -N; i < N+1; i ++)
	average = 0 // initialize
	pos = 0 // counter for the number of data points for positive velocity; initialize
	neg = 0 // counter for the number of data points for negative velocity; initialize
	name = "histResult" + num2str(i) // make an updated wave name for each different external force
	Make/o/n=1 $name=0 // set the wave to 0; this wave is used to store histogram result
	w[10] = i // set the external force
	meansOfMeansMatrix(tott, maxstep,num,N) // get the mean position
	for(z = 0; z < maxit; z++) // make 1d wave for tp 
		tpCopy[z] = tp[z][0]
	endfor
	Differentiate M /X = tpCopy// get the velocity while set the xwave to tpCopy
	for(j = 20000; j < maxit; j++)
		stableM[j - 20000] = M[j]
		if(stableM[j - 20000] > 0) // count the positve and negative velocity
			average[0] += stableM[j - 20000]
			pos++
		else
			average[1] += stableM[j - 20000]
			neg++
		endif
	endfor
	average[0] = average[0]/pos // mean value for negative velocity 
	average[1] = average[1]/neg // mean value for positive velocity 
	print(average[0])
	print(average[1])
	Histogram/B={-100,100,3} stableM,$name 
	appendToGraph average vs bin // append the average velocity to the graph (the blue and red markers)
	p = (neg/380000)*100 // the percentage of negative velocity 
	q = (pos/380000)*100 // the percentage of positve velocity 
	print(p)
	print(q)
	ModifyGraph mode(average)=3,marker(average)=19,msize(average)=5
	ModifyGraph rgb(average[0])=(16385,16388,65535)
	Display $name 
	Label bottom "Kinetechore Velocity"
	Label left "Number of data points"
	ModifyGraph mode=5
	ModifyGraph hbFill=4
	Legend/C/N=name/J/A=MC "\\s('" + name + "') left = " + num2str(p)+ "%"+"\r\\s('"+ name + "') right = " + num2str(q)+ "%"+ "\r\\s(average[0]) mean = "+ num2str(average[1]) + "\r\\s(average) mean = " + num2str(average[0])  
	print("i = " + num2str(i))
endfor

end function
//-------------------------------------------------------------------------------------------------------------------------------------------------------
function getPeakV3(tott,maxstep,num,N,target)
variable tott, maxstep, num, N
wave target
wave xp, tp
variable maxit = round(tott/maxstep) + 1
variable i, j, z
string name

//call functions to get new set of data
//meansOfMeansMatrix(tott,maxstep,num,N)
DoWindow /K peak // make the graph 
Display /N = peak

make/o/n = (maxit) tpCopy = 0
make/o/n = (maxit) velocity = 0 
make/o/n = (maxit) trace = 0 

for(i = 0; i < maxit; i ++)
	 tpCopy[i] = tp[i][0]
endfor 
trace = target // substitute wave
AppendToGraph trace vs tp[][0]
Differentiate trace /X = tpCopy /D = velocity // note: derivative divided by timestep

variable flag, currFlag, finalNum // initialize local variables
flag = 0
currFlag = 0 // counter for (wave) curr
finalNum = 0
make/o/n = (1, 3) curr = 0 // record the time, position, current distance, and flag of undecided switch
make/o/n = (10000, 3) final = nan // record the time, position, and flag of the true switch

curr[0][0] = 0.001
curr[0][1] = 0
curr[0][2] = 0
if(velocity[0] > 0)
	currFlag = 1 
else
	currFlag = -1 
endif

for(i = 1; i < maxit; i ++)
	curr[0][2] = abs(trace[i] - curr[0][1])
	if(curr[0][2] > 1) // try bigger
		if (trace[i] > curr[0][1])
			flag = 1 // assembly 
		else
			flag = -1 
		endif		
		if(finalNum == 0 || flag != final[finalNum - 1][2])
			final[finalNum][0] = curr[0][0]
			final[finalNum][1] = curr[0][1]
			final[finalNum][2] = flag
			finalNum++	
		endif 	 
	endif	
 
	if(currFlag == 1 && velocity[i] < 0) // add undecided switch to curr list if we have a switch from assembly to disassembly. 
		currFlag = -1 
		if (finalNum != 0 && currFlag != final[finalNum - 1][2])
			if (trace[i] > curr[0][1])
				curr[0][0] = tpCopy[i]
				curr[0][1] = trace[i]
				curr[0][2] = 0
			endif
		endif
	elseif (currFlag == -1 && velocity[i] > 0)  // add undecided switch to curr list if we have a switch from disassembly to assembly. 
		currFlag = 1
		if (finalNum != 0 && currFlag != final[finalNum - 1][2])
			if (trace[i] < curr[0][1])
				curr[0][0] = tpCopy[i]
				curr[0][1] = trace[i]
				curr[0][2] = 0
			endif
		endif  
	endif
	
endfor 

make/o/n = (finalNum, 2) rescue 
make/o/n = (finalNum, 2) cata 
rescue = nan
cata = nan
j = 0
z = 0
for(i = 0; i < finalNum; i ++)
	if (final[i][2] == 1)
		rescue[j][0] = final[i][0]
		rescue[j][1] = final[i][1]
		j++
	else
		cata[z][0] = final[i][0]
		cata[z][1] = final[i][1]
		z++
	endif 
endfor 
appendtoGraph rescue[][1] vs rescue [][0]
appendtoGraph cata[][1] vs cata[][0]
ModifyGraph mode(rescue)=3, marker(rescue)=1,mrkThick(rescue)=1.6,rgb(rescue)=(16385,28398,65535)
ModifyGraph mode(cata)=3 ,marker(cata)=1,mrkThick(cata)=1.6,rgb(cata)=(1,39321,19939)
Label bottom "Time[s]"
Label left "Position[um]"
Legend/C/N=text0/J "\\s(trace) Kinetochore position \r\\s(rescue) Rescue\r\\s(cata) Catastrophe"


                
end function
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// helper function
function clean(currNum)
variable currNum
wave curr
variable i

for(i = 1; i < currNum; i ++)
 curr[i][] = 0
endfor

end function
//----------------------------------------------------------------------------------------------------------------------------------------------------------------
function getPeakV2(tott, maxstep, num,N)
variable tott, maxstep, num, N
wave xp, tp, M, DM, tpCopy
variable maxit = round(tott/maxstep) + 1
variable i, j, z // counter

//call functions to get new set of data
//meansOfMeansMatrix(tott,maxstep,num,N)

DoWindow /K mean_fit // make the graph 
Display /N = mean_fit
//for(i = 0; i < N; i++)
	//appendToGraph xp[][i] vs tp[][i]
//endfor
AppendToGraph M vs tp[][0]
Duplicate/o M, DM
Differentiate DM // note: derivative divided by timestep

variable flag, true, k, currNum // initialize local variables
k = 0 // counter for (wave) final
currNum = 0 // counter for (wave) curr
make/o/n = (10000, 4) curr // record the time, position, current distance, and flag of undecided switch
make/o/n = (10000, 3) final = nan // record the time, position, and flag of the true switch
make/o/n = (10000, 1) finalT // record the time of the true switch
curr = 0
final = 0
finalT = 0

for(i = 0; i < maxit; i ++)
 tpCopy[i] = tp[i][0]
endfor

// record the first point in (wave) curr
curr[0][0] = 0.001 // time
curr[0][1] = 0 // position
curr[0][2] = 0 // current distance
if(DM[0] >= 0) // flag
	curr[0][3] = 1 // flag = 1: assembly 
	flag = 1
else 
	curr[0][3] = -1 // flag = -1: disassembly
	flag = -1
endif
currNum++

for(i = 0; i < maxit; i ++)
	for(j = 0; j < currNum; j ++)
		if(curr[j][3] == 1) // if the point should assembly 
			curr[j][2] = M[i] - curr[j][1] // update current position
		else // if the point should disassembly 
			curr[j][2] = curr[j][1] - M[i]// update current position
		endif
		
		if(curr[j][2] >= 0.1) // check for 100 nm cutoff
			FindValue /V = (curr[j][0]) finalT // check whether it is a duplicate value
			
			if (k == 0) // check whether it is true assebly/disassembly
				true = 1
			elseif (k > 0 && final[k-1][2] != curr[j][3])
				true = 1
			else
				true = 0
			endif
			
			if (true == 1 && V_value < 0) // V_value: FindValue
				final[k][0] = curr[j][0] // record time in (wave) final
				final[k][1] = curr[j][1] // record position in (wave) final 
				final[k][2] = curr[j][3] // record flag in (wave) final 
				finalT[k] = final[k][0] // record time in (wave) finalT
				curr[0][0] = final[k][0] // record the info of the true switch at the first row
				curr[0][1] = final[k][1]
				curr[0][2] = M[i] - final[k][1] // update current position
				curr[0][3] = final[k][2]
				clean(currNum) // discard all other fake switch
				currNum = 1
				k++
			endif
		endif
	endfor 
	
	if(flag == 1 && DM[i] < 0) // add undecided switch to curr list if we have a switch from assembly to disassembly. 
		curr[currNum][0] = tpCopy[i] // record time
		curr[currNum][1] = M[i] // record position
		curr[currNum][2] = 0 // record current position
		curr[currNum][3] = -1 // record flag
		flag = -1
		currNum++ // update counter
	elseif (flag == -1 && DM[i] > 0)  // add undecided switch to curr list if we have a switch from disassembly to assembly. 
		curr[currNum][0] = tpCopy[i] // record time
		curr[currNum][1] = M[i] // record position
		curr[currNum][2] = 0 // record current position
		curr[currNum][3] = 1 // record flag
		flag = 1
		currNum++
	endif
endfor

make/o/n = (k, 2) rescue 
make/o/n = (k, 2) cata 
rescue = nan
cata = nan
j = 0
z = 0
for(i = 0; i < k; i ++)
	if (final[i][2] == 1)
		rescue[j][0] = final[i][0]
		rescue[j][1] = final[i][1]
		j++
	else
		cata[z][0] = final[i][0]
		cata[z][1] = final[i][1]
		z++
	endif 
endfor 

appendtoGraph rescue[][1] vs rescue [][0]
appendtoGraph cata[][1] vs cata[][0]
ModifyGraph mode(rescue)=3, marker(rescue)=1,mrkThick(rescue)=1.6,rgb(rescue)=(16385,28398,65535)
ModifyGraph mode(cata)=3 ,marker(cata)=1,mrkThick(cata)=1.6,rgb(cata)=(1,39321,19939)
Label bottom "Time[s]"
Label left "Position[um]"
Legend/C/N=text0/J "\\s(M) Mean position\r\\s(rescue) Rescue\r\\s(cata) Catastrophe"
end function
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
// the function returns the meansOfMeanMatrix and widthMatrix that will be used in getContour function
function getMatrix(tott,maxstep,num,N)
variable tott,maxstep,num,N  // import variable for kpsim_rep function
wave o,xp,tp
variable maxit = round(tott/maxstep) + 1  // max possible steps
variable i,j,k
variable width,t, rowSum

o[1] = 0   
o[2] = 1
o[7] = 0   
o[9] = 1
o[10] = 0
make/o/n = (maxit) Mwidth // 1-D for storing width at each time point
make/o/n = (N) curr  // 1-D for storing current position of each MT
make/o/n = (maxit) M

for(i = 0; i < maxit; i ++)  // Initilize M and Mwidth to 0. If we do not initilize M and Mwidth to 0, previous value will affect the results because of line 36 and 37.
	M[i] = 0
	Mwidth[i] = 0
endfor

for (i = 0; i < num; i ++)
	kpsim_initialize(tott,maxstep,N)
	kpsim(tott,maxstep,N)
	for (j = 0; j < maxit; j ++)
		rowSum = 0
		for (k = 0; k < N; k ++)
			curr[k] = xp[j][k]         // record the current position of each tip at jth time
			rowSum = rowSum + xp[j][k]	
		endfor
		width = wavemax(curr) - wavemin(curr)  // the width of current tip 
		Mwidth[j] = (Mwidth[j] + width) // total width after num number simulations
		M[j] = (M[j] + rowSum/N)        // total mean after num number simulations
	endfor
endfor

for (i = 0; i < maxit; i ++)
	Mwidth[i] = Mwidth[i]/num   // the average width by dividing num 
	M[i] = M[i]/num             // the average mean by dividing num 
endfor

end function
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
// Call the getMatrix function and obtain the matrix for the width and means of mean contour graph. The function automatically display the width contour graph. Need to 
// manually set up the means of mean contour grpah. 
function getContour(tott,maxstep,num,N)
variable tott,maxstep,num,N // variables for calling kpsim function
variable i,j,counter // index and couter
wave w // used for adjusting external force and spring constant
variable force, spring // set the upper limit of the spring constant and external force

counter = 0
spring = 41
force = 100  // 10 simulations, higher

Make/O/D/N = (spring,force) width2D // create a matrix; row: spring constant; column: total force
make/o/n = (10000,3) mean2D
SetScale x,0,spring,width2D // spring constant: from 0 - 9 pN/um
SetScale y,0,force,width2D // total force: from 0 - 11 pN

for (i = 0; i < spring; i ++) // i+= 2
	w[8] = i
	for (j = 0; j < force; j ++)
		w[10] = j
		if (i > 15 && j < 9)    // i: spring constant; j: total force
			getMatrix(tott,.01,num,N)
		elseif (i > 30 && j > 8)
			getMatrix(tott,.01,num,N)
		else
			getMatrix(tott,maxstep,num,N) // maxstep: .1 
		endif
		wave Mwidth,M
		mean2D[counter][0] = i
		mean2D[counter][1] = j
		mean2D[counter][2] = mean(M)
		counter += 1
		width2D[i][j] = mean(Mwidth,200/maxstep,(tott/maxstep)) //calculate the mean midth after the MT reach a stable state
		print("i = " + num2str(i) + "; j = " + num2str(j)) 
	endfor 
endfor

DoWindow /K width_contour
Display /N = width_contour
AppendMatrixContour width2D

Label left "Total External Force[pN]"
Label bottom "Spring Constant, k[pN/um]"
end function
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------  
// Display the contour graph. Obtain the relationship between the bundle size and the external force for means of mean. 
function bundleContourMeans(tott,maxstep,num)
variable tott,maxstep,num // variables for calling kpsim function
variable i,j,k,counter // index and couter
wave w // used for adjusting external force and spring constant
variable force,N, forceMT

counter = 0
force = 101 // external force range: 0 - 100 pN
N = 40 // MT size range: 0 - 40
make/o/n = (4101,3) bundleMeans2D0 // 2D matrix when spring constant = 0 pN/um
make/o/n = (4101,3) bundleMeans2D40 // 2D matrix when spring constant = 40 pN/um

w[8] = 0 // set spring contant to 0 pN/um
for (i = 1; i < N; i ++)  // i = bundle size
	forceMT = i*10
	for (j = 0; j < forceMT+1; j += 3) // j = external force. Use an increment of 3. 
		w[10] = j 
		print("i = " + num2str(i) + "; j = " + num2str(j))
		meansOfMeansMatrix(tott,maxstep,num,i) // N = i
		wave M // 1-D wave: store the mean position under each time point
		bundleMeans2D0[counter][0] = i
		bundleMeans2D0[counter][1] = j/i // force per MT
		bundleMeans2D0[counter][2] = mean(M)
		counter += 1
	endfor 
endfor 

counter = 0
w[8] = 40 // set spring contant to 40 pN/um
for (i = 1; i < N; i ++)  
	forceMT = i*10 
	for (j = 0; j < forceMT+1; j += 3)  
		w[10] = j 
		print("i = " + num2str(i) + "; j = " + num2str(j))
		if (j < 50)
			meansOfMeansMatrix(tott,0.01,num,i)
		else 
			meansOfMeansMatrix(tott,maxstep,num,i)
		endif
		wave M // 1-D wave: store the mean position under each time point
		bundleMeans2D40[counter][0] = i
		bundleMeans2D40[counter][1] = j/i // force per MT
		bundleMeans2D40[counter][2] = mean(M)
		counter += 1
	endfor 
endfor 

end function
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Display the contour graph. Obtain the relathionship between the bundle size and the external force for the width. 
function bundleContourWidth(tott,maxstep,num)
variable tott,maxstep,num // variables for calling kpsim function
variable i,j,k,counter // index and couter
wave w // used for adjusting external force and spring constant
variable spring,N

spring = 41 // spring constant range: 0 - 40 pN/um 
N = 36 // MT size range: 0 - 40
make/o/n = (1435,3) bundleWidth2D0 = 0
make/o/n = (1435,3) bundleWidth2D100 = 0

counter = 0
w[10] = 0 
for (i = 1; i < N; i ++)  // i = bundle size
	for (j = 0; j < spring; j ++) // j = external force
		w[8] = j 
		print("i = " + num2str(i) + "; j = " + num2str(j))
		if (j > 13)
			widthMatrix(tott,0.01,num,i)
		else 
			widthMatrix(tott,maxstep,num,i)
		endif
		wave Mwidth 
		bundleWidth2D0[counter][0] = i
		bundleWidth2D0[counter][1] = j
		bundleWidth2D0[counter][2] = mean(Mwidth,200/maxstep,(tott/maxstep))
		counter += 1
	endfor 
endfor 

counter = 0 
for (i = 1; i < N; i ++) 
	w[10] = i*10 // Set the external force to 10 pN/MT
	for (j = 0; j < spring; j ++)  
		w[8] = j 
		print("i = " + num2str(i) + "; j = " + num2str(j))
		widthMatrix(tott,maxstep,num,i)
		wave Mwidth 
		bundleWidth2D100[counter][0] = i 
		bundleWidth2D100[counter][1] = j 
		bundleWidth2D100[counter][2] = mean(Mwidth,200/maxstep,(tott/maxstep))
		counter += 1
	endfor 
endfor 
end function
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------
// The function returns the graph of mean position and width change of the MT bundle.
function widthGraph(tott,maxstep,num,N)
variable tott,maxstep,num,N  // import variable for kpsim_rep function
wave xp,tp,o,M 
variable i

wave Mwidth,M

DoWindow /K MT_width
Display /N = MT_width
appendToGraph Mwidth vs tp[][0] // graph width vs acutal time
appendToGraph/C = (0,0,0) M vs tp[][0]
ModifyGraph lsize(M) = 5
ModifyGraph rgb(Mwidth) = (65535,0,0,30000), lsize(Mwidth) = 5

for (i = 0; i < N; i ++) // compare width with xp vs tp 
	appendToGraph/C = (24600*i, 9000*i, 12000*i), xp[][i] vs tp[][i]
endfor

String legendText = "\s(Mwidth)" + "width" + "\r\s(M)" + "mean position"

for(i = 0; i < N; i ++)
	legendText += "\r\s(xp#" + num2str(i) + ") MT" + num2str(i) 
endfor

Legend/W=MT_width legendText
ModifyGraph mirror=2
ModifyGraph fStyle=1
ModifyGraph axThick=1.5
Label left "Position(um)"
Label bottom "Time(s)"
end function

//---------------------------------------------------------------------------------------------------------------------------------------------------------------
// Obtain the 2D matrix for means of mean contour. Need to manually set up the graph. 
function meansOfMean_contour(tott,maxstep,num,N)
variable tott,maxstep,num,N  // import variable for kpsim_rep function
wave w,M
variable maxit = round(tott/maxstep) + 1  // max possible steps
variable i,j,counter
variable force,spring
 

make/o/n = (1000,3) mean2D
counter = 0
force = 15
spring = 20

for (i = 0; i < spring; i ++)
	w[8] = i
	for(j = 0; j < force; j ++)
		w[10] = j
		meansOfMeansMatrix(tott,maxstep,num,N)
		wave M
		mean2D[counter][0] = i
		mean2D[counter][1] = j
		mean2D[counter][2] = mean(M)
		counter += 1
	endfor
endfor
end function

//---------------------------------------------------------------------------------------------------------------------------------------------------------
// Use a the 1-D wave to store the mean position at each time point. 
function meansOfMeansMatrix(tott,maxstep,num,N)
variable tott,maxstep,num,N  // import variable for kpsim_rep function
wave xp,tp,o,M 
variable maxit = round(tott/maxstep) + 1  // max possible steps
variable i,j,k
variable rowSum 

o[1] = 0  // no single state
o[2] = 1
o[7] = 0  // no detachment
o[9] = 1
o[10] = 0

make/o/n = (maxit) M // the column of the wave records mean position at each time step under a specific forc

for(i = 0; i < maxit; i ++)  // Initilize M and Mwidth to 0. If we do not initilize M and Mwidth to 0, previous value will affect the results because of line 36 and 37.
	M[i] = 0
endfor

for (i = 0; i < num; i ++)
	kpsim_initialize(tott,maxstep,N)
	kpsim(tott,maxstep,N)
	for (j = 0; j < maxit; j ++)
		rowSum = 0
		for (k = 0; k < N; k ++)
			rowSum = rowSum + xp[j][k]	
		endfor
		M[j] = M[j] + rowSum/N        // total mean after num number simulations
	endfor
endfor

for (i = 0; i < maxit; i ++)
	M[i] = M[i]/num             // the average mean by dividing num 
endfor

end function
//----------------------------------------------------------------------------------------------------------------------------------------------------------
// Use a 1-D wave to store the width at each time point. 
function widthMatrix(tott,maxstep,num,N)
variable tott,maxstep,num,N  // import variable for kpsim_rep function
wave o,xp,tp
variable maxit = round(tott/maxstep) + 1  // max possible steps
variable i,j,k
variable width,t

o[1] = 0   
o[2] = 1
o[7] = 0   
o[9] = 1
o[10] = 0
make/o/n = (maxit) Mwidth // 1-D for storing width at each time point
make/o/n = (N) curr  // 1-D for storing current position of each MT

for (i = 0; i < maxit; i ++)
	Mwidth[i] = 0
endfor

for (i = 0; i < num; i ++)
	kpsim_initialize(tott,maxstep,N)
	kpsim(tott,maxstep,N)
	for (j = 0; j < maxit; j ++)
		for (k = 0; k < N; k ++)
			curr[k] = xp[j][k]         // record the current position of each tip at jth time
		endfor
		width = wavemax(curr) - wavemin(curr)  // the width of current tip 
		Mwidth[j] = Mwidth[j] + width // total width after num number simulations
	endfor
endfor

for (i = 0; i < maxit; i ++)
	Mwidth[i] = Mwidth[i]/num   // the average width by dividing num 
endfor

end function