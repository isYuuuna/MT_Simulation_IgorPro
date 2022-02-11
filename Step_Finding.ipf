#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

function fullpackage()
wave o

o[8] = 8// kappa
o[10] = 10 // total force
kpsim_rep(1000,0.01,1,5)
getKTvelocity(1000,0.01,5)
wave KTp, KTv
GaussianFilter(1000,0.01)
wave filtered
decimate(filtered,100)
wave opDecimate
GetChangePointIndexes(opDecimate,1)
wave changePointIndexes
getPosition(changePointIndexes)
end function

//----------------------------------------------------------
function GetPosition(data)
wave data
variable i

wave opDecimate
wavestats data
make/o/n = (V_npnts) position 
for (i = 0; i < V_npnts; i ++)
	position[i] = opDecimate[data[i]]
endfor

end function

//--------------------------------------------------------------
function GetChangePointIndexes(data, minDistance)
wave data
variable minDistance

variable n = numpnts(data)
//variable penalty = 3*ln(n)
//variable penalty = 40
variable penalty = ln(n)*40
variable k = min(n,ceil(4*ln(n)))
//variable k = min(n,40)

GetPartialSums(data,k)
wave partialSums

make/o/n = (n+1) bestCost = nan
bestCost[0] = -penalty
variable currentTau
for(currentTau = minDistance; currentTau < 2*minDistance; currentTau++)
	bestCost[currentTau] = getSegmentCost(partialSums,0, currentTau,k,n)
endfor

make/o/n = (n+1) previousChangePointIndex = 0
make/o/n = (n+1) previousTaus = nan
previousTaus[0] = 0
previousTaus[1] = minDistance
make/o/n = (n+1) costForPreviousTau = nan

variable num = 0

for(currentTau = 2 * minDistance; currentTau < n + 1; currentTau++)
	costForPreviousTau = nan
	variable costForPreviousTauCount = 0
	
	variable previousTau,i
	wavestats previousTaus
	for(i = 0; i < V_npnts; i ++)
		previousTau = previousTaus[i]
		num = getSegmentCost(partialSums,previousTau,currentTau,k,n)
		costForPreviousTau[costForPreviousTauCount] = bestCost[previousTau] + num + penalty
		costForPreviousTauCount++
	endfor
	
	variable bestPreviousTauIndex = WhichMin(costForPreviousTau)
	bestCost[currentTau] = costForPreviousTau[bestPreviousTauIndex]
    previousChangePointIndex[currentTau] = previousTaus[bestPreviousTauIndex]
    
    variable currentBestCost = bestCost[currentTau]
    variable newPreviousTausSize = 0
    variable j    
    waveStats previousTaus    
    for(i = 0; i < V_npnts; i ++)
    	if (costForPreviousTau[i] < currentBestCost + penalty)
        	previousTaus[newPreviousTausSize] = previousTaus[i]
        	newPreviousTausSize++
        endif
    endfor
    
    waveStats previousTaus 
    variable removeCount = V_npnts - newPreviousTausSize
    removeRange(previousTaus, newPreviousTausSize, removeCount, n)
    wave previousTaus
    
     waveStats previousTaus
     previousTaus[V_npnts] = currentTau - (minDistance - 1)   
endfor

make/o/n = (n) changePointIndexes = nan
variable currentIndex = previousChangePointIndex[n]
do
	waveStats changePointIndexes
	changePointIndexes[V_npnts] = currentIndex - 1
	currentIndex = previousChangePointIndex[currentIndex]
while(currentIndex != 0)

sort changePointIndexes,changePointIndexes

end function

//----------------------------------------------------
function getPartialSums(data, k)
wave data
variable k

variable n = numpnts(data)
make/o/n = (k,(n + 1)) partialSums = 0
Duplicate/O data, sortedData
sort sortedData,sortedData

variable i
for(i = 0; i < k; i++)
	variable z = -1 + (2 * i + 1)/k
	variable p = 1 / (1 + (2 * n - 1)^(-z))
	variable t = sortedData[Trunc((n-1)*p)]
	
	variable tau
	for(tau = 1; tau <= n; tau++)
		partialSums[i][tau] = partialSums[i][tau - 1];
		if (data[tau - 1] < t)
			partialSums[i][tau] += 2
		endif
		if (data[tau - 1] == t)
			partialSums[i][tau] += 1
		endif
	endfor	
endfor 

end function

//---------------------------------------------------------
function GetSegmentCost(partialSums, tau1, tau2, k, n)
wave partialSums
variable tau1, tau2, k, n

variable summ = 0
variable i
for(i = 0; i < k; i ++)
	variable actualSum = partialSums[i][tau2] - partialSums[i][tau1]
	
	if (actualSum != 0 && actualSum != (tau2 - tau1)*2)
		variable fit = actualSum * 0.5/ (tau2 - tau1)
		variable lnp = (tau2 - tau1)*(fit * ln(fit) + (1 - fit)*ln(1 - fit))
		summ += lnp
	endif
	
endfor

variable c = -ln(2*n - 1)
return 2*c / k*summ
end function 
//------------------------------------------------------------------
function WhichMin(values)
wave values

variable minValue = values[0]
variable minIndex = 0
variable i

wavestats values
for (i = 1; i < V_npnts; i ++)
	if(values[i] < minValue)
		minvalue = values[i]
		minIndex = i
	endif
endfor

return minIndex
end function

//-------------------------------------------------------------------
function removeRange(previousTaus, start, count, n)
wave previousTaus
variable start, count, n
variable i

make/o/n = (n) previousTausCopy = previousTaus

wavestats previousTaus
variable length = V_npnts
for (i = start; i < count; i ++)
	previousTaus[i] = previousTausCopy[i + count]
endfor

if(length != 0)
	for (i = length - count; i < length; i ++)
		previousTaus[i] = nan
	endfor  
endif 

end function